# -*- coding: utf-8 -*-
# YAFF is yet another force-field code
# Copyright (C) 2011 - 2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
# Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>, Center for Molecular Modeling
# (CMM), Ghent University, Ghent, Belgium; all rights reserved unless otherwise
# stated.
#
# This file is part of YAFF.
#
# YAFF is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# YAFF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--
'''Barostats'''


import numpy as np

from molmod import boltzmann, femtosecond, kjmol, bar, atm

from yaff.log import log, timer
from yaff.sampling.utils import get_random_vel, cell_symmetrize, get_random_vel_press, \
    get_ndof_internal_md, clean_momenta, get_ndof_baro
from yaff.sampling.verlet import VerletHook
from yaff.sampling.iterative import StateItem


__all__ = [
    'TBCombination', 'McDonaldBarostat', 'BerendsenBarostat', 'LangevinBarostat',
    'MTKBarostat'
]

class TBCombination(VerletHook):
    def __init__(self, thermostat, barostat, start = 0):
        """
            VerletHook combining an arbitrary Thermostat and Barostat instance, which
            ensures these instances are called in the correct succession, and possible
            coupling between both is handled correctly.

            **Arguments:**

            thermostat
                A Thermostat instance

            barostat
                A Barostat instance

        """
        self.thermostat = thermostat
        self.barostat = barostat
        self.start = start
        # verify if thermostat and barostat instances are currently supported in yaff
        if not self.verify():
            self.barostat = thermostat
            self.thermostat = barostat
            if not self.verify():
                raise TypeError('The Thermostat or Barostat instance is not supported (yet)')
        self.step_thermo = self.thermostat.step
        self.step_baro = self.barostat.step
        VerletHook.__init__(self, start, min(self.step_thermo, self.step_baro))

    def init(self, iterative):
        # initialize the thermostat and barostat separately
        self.thermostat.init(iterative)
        self.barostat.init(iterative)
        # variables which will determine the coupling between thermostat and barostat
        self.chainvel0 = None
        self.G1_add = None

    def pre(self, iterative):
        # determine whether the barostat should be called
        if self.expectscall(iterative, 'baro'):
            from yaff.sampling.nvt import NHCThermostat
            if isinstance(self.thermostat, NHCThermostat):
                # in case the barostat is coupled with a NHC thermostat:
                # v_{xi,1} is needed to update v_g
                self.chainvel0 = self.thermostat.chain.vel[0]
            # actual barostat update
            self.barostat.pre(iterative, self.chainvel0)
        # determine whether the thermostat should be called
        if self.expectscall(iterative, 'thermo'):
            if isinstance(self.barostat, MTKBarostat):
                # in case the thermostat is coupled with a MTK barostat:
                # update equation of v_{xi,1} is altered via G_1
                self.G1_add = self.barostat.add_press_cont()
            # actual thermostat update
            self.thermostat.pre(iterative, self.G1_add)

    def post(self, iterative):
        # determine whether the thermostat should be called
        if self.expectscall(iterative, 'thermo'):
            if isinstance(self.barostat, MTKBarostat):
                # in case the thermostat is coupled with a MTK barostat:
                # update equation of v_{xi,1} is altered via G_1
                self.G1_add = self.barostat.add_press_cont()
            # actual thermostat update
            self.thermostat.post(iterative, self.G1_add)
        # determine whether the barostat should be called
        if self.expectscall(iterative, 'baro'):
            from yaff.sampling.nvt import NHCThermostat, LangevinThermostat
            if isinstance(self.thermostat, NHCThermostat):
                # in case the barostat is coupled with a NHC thermostat:
                # v_{xi,1} is needed to update v_g
                self.chainvel0 = self.thermostat.chain.vel[0]
            #if isinstance(self.thermostat, LangevinThermostat):
                # Langevin thermostat: correct for final ekin and epot (1)
                #self.econs_correction += iterative.ekin + iterative.epot
            # actual barostat update
            self.barostat.post(iterative, self.chainvel0)
           # if isinstance(self.thermostat, LangevinThermostat):
                # Langevin thermostat: correct for final ekin and epot (2)
                #self.econs_correction += -iterative.ekin - iterative.epot
        # update the correction on E_cons due to thermostat and barostat
        self.econs_correction = self.thermostat.econs_correction + self.barostat.econs_correction
        if isinstance(self.thermostat, NHCThermostat):
            if isinstance(self.barostat, MTKBarostat) and self.barostat.baro_thermo is not None: pass
            else:
                # extra correction necessary if particle NHC thermostat is used to thermostat the barostat
                kt = boltzmann*self.thermostat.temp
                baro_ndof = self.barostat.baro_ndof
                self.econs_correction += baro_ndof*kt*self.thermostat.chain.pos[0]

    def expectscall(self, iterative, kind):
        # returns whether the thermostat/barostat should be called in this iteration
        if kind == 'thermo':
            return iterative.counter >= self.start and (iterative.counter - self.start) % self.step_thermo == 0
        if kind == 'baro':
            return iterative.counter >= self.start and (iterative.counter - self.start) % self.step_baro == 0

    def verify(self):
        # returns whether the thermostat and barostat instances are currently supported by yaff
        from yaff.sampling.nvt import AndersenThermostat, NHCThermostat, LangevinThermostat, BerendsenThermostat
        thermo_correct = False
        baro_correct = False
        thermo_list = [AndersenThermostat, NHCThermostat, LangevinThermostat, BerendsenThermostat]
        baro_list = [McDonaldBarostat, BerendsenBarostat, LangevinBarostat, MTKBarostat]
        if any(isinstance(self.thermostat, thermo) for thermo in thermo_list):
            thermo_correct = True
        if any(isinstance(self.barostat, baro) for baro in baro_list):
            baro_correct = True
        return (thermo_correct and baro_correct)


class McDonaldBarostat(VerletHook):
    def __init__(self, temp, press, start=0, step=1, amp=1e-3):
        """
           Warning: this code is not fully tested yet!

           **Arguments:**

           temp
                The average temperature of the NpT ensemble

           press
                The external pressure of the NpT ensemble

           **Optional arguments:**

           start
                The first iteration at which this hook is called

           step
                The number of iterations between two subsequent calls to this
                hook.

           amp
                The amplitude of the changes in the logarithm of the volume.
        """
        self.temp = temp
        self.press = press
        self.amp = amp
        self.dim = 3
        self.baro_ndof = 1
        VerletHook.__init__(self, start, step)

    def init(self, iterative):
        pass

    def pre(self, iterative, chainvel0 = None):
        pass

    def post(self, iterative, chainvel0 = None):
        def compute(pos, rvecs):
            iterative.pos[:] = pos
            iterative.gpos[:] = 0.0
            iterative.ff.update_rvecs(rvecs)
            iterative.ff.update_pos(pos)
            iterative.epot = iterative.ff.compute(iterative.gpos)
            iterative.acc = -iterative.gpos/iterative.masses.reshape(-1,1)

        natom = iterative.ff.system.natom
        with timer.section('AMB'):
            # A) Change the logarithm of the volume isotropically.
            scale = np.exp(np.random.uniform(-self.amp, self.amp))
            # A.0) Keep track of old state
            vol0 = iterative.ff.system.cell.volume
            epot0 = iterative.epot
            rvecs0 = iterative.ff.system.cell.rvecs.copy()
            pos0 = iterative.pos.copy()
            # A.1) scale the system and recompute the energy
            compute(pos0*scale, rvecs0*scale)
            epot1 = iterative.epot
            vol1 = iterative.ff.system.cell.volume
            # A.2) compute the acceptance ratio
            beta = 1/(boltzmann*self.temp)
            arg = (epot1 - epot0 + self.press*(vol1 - vol0) - (natom+1)/beta*np.log(vol1/vol0))
            accepted = arg < 0 or np.random.uniform(0, 1) < np.exp(-beta*arg)
            if accepted:
                # add a correction to the conserved quantity
                self.econs_correction += epot0 - epot1
            else:
                # revert the cell and the positions in the original state
                compute(pos0, rvecs0)
            # B) Change the velocities
            ekin0 = iterative._compute_ekin()
            iterative.vel[:] = get_random_vel(self.temp, False, iterative.masses)
            # C) Update the kinetic energy and the reference for the conserved quantity
            ekin1 = iterative._compute_ekin()
            self.econs_correction += ekin0 - ekin1
            if log.do_medium:
                with log.section('AMB'):
                    s = {True: 'accepted', False: 'rejected'}[accepted]
                    log('BARO   volscale %10.7f      arg %s      %s' % (scale, log.energy(arg), s))
                    if accepted:
                        log('BARO   energy change %s      (new vol)**(1/3) %s' % (
                            log.energy(epot1 - epot0), log.length(vol1**(1.0/3.0))
                        ))
                    log('THERMO energy change %s' % log.energy(ekin0 - ekin1))



class BerendsenBarostat(VerletHook):
    def __init__(self, ff, temp, press, start=0, step=1, timecon=1000*femtosecond, beta=4.57e-5/bar, anisotropic=True, vol_constraint=False):
        """
            This hook implements the Berendsen barostat. The equations are derived in:

                Berendsen, H. J. C.; Postma, J. P. M.; van Gunsteren, W. F.;
                Dinola, A.; Haak, J. R. J. Chem. Phys. 1984, 81, 3684-3690

            **Arguments:**

            ff
                A ForceField instance.

            temp
                The temperature of thermostat.

            press
                The applied pressure for the barostat.

            **Optional arguments:**

            start
                The step at which the thermostat becomes active.

            timecon
                The time constant of the Berendsen barostat.

            beta
                The isothermal compressibility, conventionally the compressibility of liquid water

            anisotropic
                Defines whether anisotropic cell fluctuations are allowed.

            vol_constraint
                Defines whether the volume is allowed to fluctuate.
        """
        self.temp = temp
        self.press = press
        self.mass_press = 3.0*timecon/beta
        self.anisotropic = anisotropic
        self.vol_constraint = vol_constraint
        self.dim = ff.system.cell.nvec
        # determine the number of degrees of freedom associated with the unit cell
        self.baro_ndof = get_ndof_baro(self.dim, self.anisotropic, self.vol_constraint)
        if self.anisotropic:
            # symmetrize the cell tensor
            cell_symmetrize(ff)
        self.cell = ff.system.cell.rvecs.copy()
        self.dim = ff.system.cell.nvec
        VerletHook.__init__(self, start, 1)

    def init(self, iterative):
        self.timestep_press = iterative.timestep
        # compute gpos and vtens, since they differ
        # after symmetrising the cell tensor
        iterative.gpos[:] = 0.0
        iterative.vtens[:] = 0.0
        iterative.epot = iterative.ff.compute(iterative.gpos,iterative.vtens)
        if iterative.ndof is None:
            iterative.ndof = get_ndof_internal_md(iterative.pos.shape[0], iterative.ff.system.cell.nvec)
        # rescaling of the barostat mass, to be in accordance with Langevin and MTK
        self.mass_press *= np.sqrt(iterative.ndof)

    def pre(self, iterative, chainvel0 = None):
        pass

    def post(self, iterative, chainvel0 = None):
        # for bookkeeping purposes
        epot0 = iterative.epot
        # calculation of the internal pressure tensor
        ptens = (np.dot(iterative.vel.T*iterative.masses, iterative.vel) - iterative.vtens)/iterative.ff.system.cell.volume
        # determination of mu
        dmu = self.timestep_press/self.mass_press*(self.press*np.eye(3)-ptens)
        if self.vol_constraint:
            dmu -= np.trace(dmu)/self.dim*np.eye(self.dim)
        mu = np.eye(3) - dmu
        mu = 0.5*(mu+mu.T)
        if not self.anisotropic:
            mu = ((np.trace(mu)/3.0)**(1.0/3.0))*np.eye(3)
        # updating the positions and cell vectors
        pos_new = np.dot(iterative.pos, mu)
        rvecs_new = np.dot(iterative.rvecs, mu)
        iterative.ff.update_pos(pos_new)
        iterative.pos[:] = pos_new
        iterative.ff.update_rvecs(rvecs_new)
        iterative.rvecs[:] = rvecs_new
        # calculation of the virial tensor
        iterative.gpos[:] = 0.0
        iterative.vtens[:] = 0.0
        iterative.epot = iterative.ff.compute(iterative.gpos,iterative.vtens)
        epot1 = iterative.epot
        self.econs_correction += epot0 - epot1


class LangevinBarostat(VerletHook):
    def __init__(self, ff, temp, press = 1e6, start=0, step=1, timecon=1000*femtosecond, anisotropic=True, vol_constraint=False):
        """
            This hook implements the Langevin barostat. The equations are derived in:

                Feller, S. E.; Zhang, Y.; Pastor, R. W.; Brooks, B. R.
                J. Chem. Phys. 1995, 103, 4613-4621

            **Arguments:**

            ff
                A ForceField instance.

            temp
                The temperature of thermostat.

            press
                The applied pressure for the barostat.

            **Optional arguments:**

            start
                The step at which the barostat becomes active.

            timecon
                The time constant of the Langevin barostat.

            anisotropic
                Defines whether anisotropic cell fluctuations are allowed.

            vol_constraint
                Defines whether the volume is allowed to fluctuate.
        """
        self.temp = temp
        self.press = press
        self.timecon = timecon
        self.anisotropic = anisotropic
        self.vol_constraint = vol_constraint
        self.dim = ff.system.cell.nvec
        # determine the number of degrees of freedom associated with the unit cell
        self.baro_ndof = get_ndof_baro(self.dim, self.anisotropic, self.vol_constraint)
        if self.anisotropic:
            # symmetrize the cell tensor
            cell_symmetrize(ff)
        self.cell = ff.system.cell.rvecs.copy()
        VerletHook.__init__(self, start, step)

    def init(self, iterative):
        self.timestep_press = iterative.timestep
        clean_momenta(iterative.pos, iterative.vel, iterative.masses, iterative.ff.system.cell)
        # define the barostat 'mass'
        self.ndof = get_ndof_internal_md(len(iterative.ff.system.numbers), iterative.ff.system.cell.nvec)
        self.mass_press = (self.ndof+3)/3*boltzmann*self.temp*(self.timecon/(2*np.pi))**2
        # define initial barostat velocity
        self.vel_press = get_random_vel_press(self.mass_press, self.temp)
        # make sure the volume of the cell will not change if applicable
        if self.vol_constraint:
            self.vel_press -= np.trace(self.vel_press)/3*np.eye(3)
        if not self.anisotropic:
            self.vel_press = self.vel_press[0][0]

        # compute gpos and vtens, since they differ
        # after symmetrising the cell tensor
        iterative.gpos[:] = 0.0
        iterative.vtens[:] = 0.0
        iterative.epot = iterative.ff.compute(iterative.gpos,iterative.vtens)

    def pre(self, iterative, chainvel0 = None):
        # bookkeeping
        epot0 = iterative.epot
        ekin0 = iterative.ekin
        # the actual update
        self.baro(iterative, chainvel0)
        # some more bookkeeping
        epot1 = iterative.epot
        ekin1 = iterative.ekin
        self.econs_correction += epot0 - epot1 + ekin0 - ekin1

    def post(self, iterative, chainvel0 = None):
        # bookkeeping
        epot0 = iterative.epot
        ekin0 = iterative.ekin
        # the actual update
        self.baro(iterative, chainvel0)
        # some more bookkeeping
        epot1 = iterative.epot
        ekin1 = iterative.ekin
        self.econs_correction += epot0 - epot1 + ekin0 - ekin1

    def baro(self, iterative, chainvel0):
        def update_baro_vel():
            # updates the barostat velocity tensor
            # iL h/(8*tau)
            self.vel_press *= np.exp(-self.timestep_press/(8*self.timecon))
            if chainvel0 is not None:
                # iL v_{xi} v_g h/8: extra contribution due to NHC thermostat
                self.vel_press *= np.exp(-self.timestep_press*chainvel0/8)
            # definition of P_intV and G
            ptens_vol = np.dot(iterative.vel.T*iterative.masses, iterative.vel) - iterative.vtens
            ptens_vol = 0.5*(ptens_vol.T + ptens_vol)
            G = (ptens_vol+(2.0*iterative.ekin/self.ndof-self.press*iterative.ff.system.cell.volume)*np.eye(3))/self.mass_press
            R = self.getR()
            if self.vol_constraint:
                G -= np.trace(G)/self.dim*np.eye(self.dim)
                R -= np.trace(R)/self.dim*np.eye(self.dim)
            if not self.anisotropic:
                G = np.trace(G)
                R = R[0][0]
            # iL (G_g-R_p/W) h/4
            self.vel_press += (G-R/self.mass_press)*self.timestep_press/4
            # iL h/(8*tau)
            self.vel_press *= np.exp(-self.timestep_press/(8*self.timecon))
            if chainvel0 is not None:
                # iL v_{xi} v_g h/8: extra contribution due to NHC thermostat
                self.vel_press *= np.exp(-self.timestep_press*chainvel0/8)

        # first part of the barostat velocity tensor update
        update_baro_vel()

        # iL v_g h/2
        if self.anisotropic:
            Dr, Qg = np.linalg.eigh(self.vel_press)
            Daccr = np.diagflat(np.exp(Dr*self.timestep_press/2))
            rot_mat = np.dot(np.dot(Qg, Daccr), Qg.T)
            pos_new = np.dot(iterative.pos, rot_mat)
            rvecs_new = np.dot(iterative.rvecs, rot_mat)
        else:
            c = np.exp(self.vel_press*self.timestep_press/2)
            pos_new = c*iterative.pos
            rvecs_new = c*iterative.rvecs

        # update the positions and cell vectors
        iterative.ff.update_pos(pos_new)
        iterative.pos[:] = pos_new
        iterative.ff.update_rvecs(rvecs_new)
        iterative.rvecs[:] = rvecs_new

        # update the potential energy
        iterative.gpos[:] = 0.0
        iterative.vtens[:] = 0.0
        iterative.epot = iterative.ff.compute(iterative.gpos,iterative.vtens)

        # -iL (v_g + Tr(v_g)/ndof) h/2
        if self.anisotropic:
            if self.vol_constraint:
                Dg, Eg = np.linalg.eigh(self.vel_press)
            else:
                Dg, Eg = np.linalg.eigh(self.vel_press+(np.trace(self.vel_press)/self.ndof)*np.eye(3))
            Daccg = np.diagflat(np.exp(-Dg*self.timestep_press/2))
            rot_mat = np.dot(np.dot(Eg, Daccg), Eg.T)
            vel_new = np.dot(iterative.vel, rot_mat)
        else:
            vel_new = np.exp(-((1.0+3.0/self.ndof)*self.vel_press)*self.timestep_press/2) * iterative.vel
        iterative.vel[:] = vel_new

        # update kinetic energy
        iterative.ekin = iterative._compute_ekin()

        # second part of the barostat velocity tensor update
        update_baro_vel()

    def getR(self):
        shape = 3, 3
        # generate random 3x3 tensor
        rand = np.random.normal(0, 1, shape)*np.sqrt(2*self.mass_press*boltzmann*self.temp/(self.timestep_press*self.timecon))
        R = np.zeros(shape)
        # create initial symmetric pressure velocity tensor
        for i in xrange(3):
            for j in xrange(3):
                if i >= j:
                    R[i,j] = rand[i,j]
                else:
                    R[i,j] = rand[j,i]
        return R


class MTKBarostat(VerletHook):
    def __init__(self, ff, temp, press=1*atm, start=0, step=1, timecon=1000*femtosecond, anisotropic=True, vol_constraint=False, baro_thermo=None):
        """
            This hook implements the combination of the Nosé-Hoover chain thermostat
            and the Martyna-Tobias-Klein barostat. The equations are derived in:

                Martyna, G. J.; Tobias, D. J.: Klein, M. L. J. Chem. Phys. 1994,
                101, 4177-4189.

            The implementation (used here) of a symplectic integrator of this thermostat
            and barostat is discussed in

                Martyna, G. J.;  Tuckerman, M. E.;  Tobias, D. J.;  Klein,
                M. L. Mol. Phys. 1996, 87, 1117-1157.

            **Arguments:**

            ff
                A ForceField instance.

            temp
                The temperature of thermostat.

            press
                The applied pressure for the barostat.

            **Optional arguments:**

            start
                The step at which the thermostat becomes active.

            timecon
                The time constant of the Martyna-Tobias-Klein barostat.

            anisotropic
                Defines whether anisotropic cell fluctuations are allowed.

            vol_constraint
                Defines whether the volume is allowed to fluctuate.

            baro_thermo
                NHCThermostat instance, coupled directly to the barostat
        """
        self.temp = temp
        self.press = press
        self.timecon_press = timecon
        self.anisotropic = anisotropic
        self.vol_constraint = vol_constraint
        self.baro_thermo = baro_thermo
        self.dim = ff.system.cell.nvec
        # determine the number of degrees of freedom associated with the unit cell
        self.baro_ndof = get_ndof_baro(self.dim, self.anisotropic, self.vol_constraint)
        if self.anisotropic:
            # symmetrize the cell tensor
            cell_symmetrize(ff)
        self.cell = ff.system.cell.rvecs.copy()
        VerletHook.__init__(self, start, step)

    def set_ndof(self, ndof):
        # allocate degrees of freedom
        angfreq = 2*np.pi/self.timecon_press
        self.mass_press = (ndof+self.dim**2)*boltzmann*self.temp/angfreq**2
        # define initial barostat velocity
        self.vel_press = get_random_vel_press(self.mass_press, self.temp)
        if not self.anisotropic:
            self.vel_press = self.vel_press[0][0]
        # initialize the barostat thermostat if present
        if self.baro_thermo is not None:
            self.baro_thermo.chain.timestep = iterative.timestep
            self.baro_thermo.chain.set_ndof(self.baro_ndof)
        # make sure the volume of the cell will not change if applicable
        if self.vol_constraint:
            self.vel_press -= np.trace(self.vel_press)/3*np.eye(3)
        # compute gpos and vtens, since they differ
        # after symmetrising the cell tensor
        iterative.gpos[:] = 0.0
        iterative.vtens[:] = 0.0
        iterative.epot = iterative.ff.compute(iterative.gpos,iterative.vtens)

    def pre(self, iterative, chainvel0 = None):
        if self.baro_thermo is not None:
            # overrule the chainvel0 argument in the TBC instance
            chainvel0 = self.baro_thermo.chain.vel[0]
        # propagate the barostat
        self.baro(iterative, chainvel0)
        # propagate the barostat thermostat if present
        if self.baro_thermo is not None:
            # determine the barostat kinetic energy
            ekin_baro = self._compute_ekin_baro()
            # update the barostat thermostat, without updating v_g (done in self.baro)
            vel_press_copy = np.zeros(self.vel_press.shape)
            vel_press_copy[:] = self.vel_press
            vel_baro, ekin_baro = self.baro_thermo.chain(ekin_baro, vel_press_copy, 0)

    def post(self, iterative, chainvel0 = None):
        # propagate the barostat thermostat if present
        if self.baro_thermo is not None:
            # determine the barostat kinetic energy
            ekin_baro = self._compute_ekin_baro()
            # update the barostat thermostat, without updating v_g (done in self.baro)
            vel_press_copy = np.zeros(self.vel_press.shape)
            vel_press_copy[:] = self.vel_press
            vel_baro, ekin_baro = self.baro_thermo.chain(ekin_baro, vel_press_copy, 0)
            # overrule the chainvel0 argument in the TBC instance
            chainvel0 = self.baro_thermo.chain.vel[0]
        # propagate the barostat
        self.baro(iterative, chainvel0)
        # calculate the correction due to the barostat alone
        self.econs_correction = self.press*iterative.ff.system.cell.volume + self._compute_ekin_baro()
        if self.baro_thermo is not None:
            # add the correction due to the barostat thermostat
            self.econs_correction += self.baro_thermo.chain.get_econs_correction()

    def baro(self, iterative, chainvel0):
        def update_baro_vel():
            # updates the barostat velocity tensor
            if chainvel0 is not None:
                # iL v_{xi} v_g h/8
                self.vel_press *= np.exp(-self.timestep_press*chainvel0/8)
            # definition of P_intV and G
            ptens_vol = np.dot(iterative.vel.T*iterative.masses, iterative.vel) - iterative.vtens
            ptens_vol = 0.5*(ptens_vol.T + ptens_vol)
            G = (ptens_vol+(2.0*iterative.ekin/self.ndof-self.press*iterative.ff.system.cell.volume)*np.eye(3))/self.mass_press
            if not self.anisotropic:
                G = np.trace(G)
            if self.vol_constraint:
                G -= np.trace(G)/self.dim*np.eye(self.dim)
            # iL G_g h/4
            self.vel_press += G*self.timestep_press/4
            if chainvel0 is not None:
                # iL v_{xi} v_g h/8
                self.vel_press *= np.exp(-self.timestep_press*chainvel0/8)

        # first part of the barostat velocity tensor update
        update_baro_vel()

        # iL v_g h/2
        if self.anisotropic:
            Dr, Qg = np.linalg.eigh(self.vel_press)
            Daccr = np.diagflat(np.exp(Dr*self.timestep_press/2))
            rot_mat = np.dot(np.dot(Qg, Daccr), Qg.T)
            pos_new = np.dot(iterative.pos, rot_mat)
            rvecs_new = np.dot(iterative.rvecs, rot_mat)
        else:
            c = np.exp(self.vel_press*self.timestep_press/2)
            pos_new = c*iterative.pos
            rvecs_new = c*iterative.rvecs

        # update the positions and cell vectors
        iterative.ff.update_pos(pos_new)
        iterative.pos[:] = pos_new
        iterative.ff.update_rvecs(rvecs_new)
        iterative.rvecs[:] = rvecs_new

        # update the potential energy
        iterative.gpos[:] = 0.0
        iterative.vtens[:] = 0.0
        iterative.epot = iterative.ff.compute(iterative.gpos,iterative.vtens)

        # -iL (v_g + Tr(v_g)/ndof) h/2
        if self.anisotropic:
            if self.vol_constraint:
                Dg, Eg = np.linalg.eigh(self.vel_press)
            else:
                Dg, Eg = np.linalg.eigh(self.vel_press+(np.trace(self.vel_press)/self.ndof)*np.eye(3))
            Daccg = np.diagflat(np.exp(-Dg*self.timestep_press/2))
            rot_mat = np.dot(np.dot(Eg, Daccg), Eg.T)
            vel_new = np.dot(iterative.vel, rot_mat)
        else:
            vel_new = np.exp(-((1.0+3.0/self.ndof)*self.vel_press)*self.timestep_press/2) * iterative.vel

        # update the velocities and the kinetic energy
        iterative.vel[:] = vel_new
        iterative.ekin = iterative._compute_ekin()

        # second part of the barostat velocity tensor update
        update_baro_vel()

    def add_press_cont(self):
        kt = self.temp*boltzmann
        # pressure contribution to thermostat: kinetic cell tensor energy
        # and extra degrees of freedom due to cell tensor
        if self.baro_thermo is None:
            return 2*self._compute_ekin_baro() - self.baro_ndof*kt
        else:
            # if a barostat thermostat is present, the thermostat is decoupled
            return 0

    def _compute_ekin_baro(self):
        # returns the kinetic energy associated with the cell fluctuations
        if self.anisotropic:
            return 0.5*self.mass_press*np.trace(np.dot(self.vel_press.T,self.vel_press))
        else:
            return 0.5*self.mass_press*self.vel_press**2
