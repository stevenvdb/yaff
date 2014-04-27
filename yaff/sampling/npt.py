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

from molmod import boltzmann, femtosecond, kjmol, bar

from yaff.log import log, timer
from yaff.sampling.utils import get_random_vel, cell_symmetrize, get_random_vel_press, \
    get_ndof_internal_md, clean_momenta
from yaff.sampling.verlet import VerletHook


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
                raise TypeError('The Thermostat or Barostat instance is not supported')
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
            self.barostat.pre(iterative, self.chainvel0)
        # determine whether the thermostat should be called
        if self.expectscall(iterative, 'thermo'):
            if isinstance(self.barostat, MTKBarostat):
                # in case the thermostat is coupled with a MTK barostat:
                # update equation of v_{xi,1} is altered via G_1
                self.G1_add = self.barostat.add_press_cont()
            self.thermostat.pre(iterative, self.G1_add)

    def post(self, iterative):
        # determine whether the thermostat should be called
        if self.expectscall(iterative, 'thermo'):
            if isinstance(self.barostat, MTKBarostat):
                # in case the thermostat is coupled with a MTK barostat:
                # update equation of v_{xi,1} is altered via G_1
                self.G1_add = self.barostat.add_press_cont()
            self.thermostat.post(iterative, self.G1_add)
        # determine whether the barostat should be called
        if self.expectscall(iterative, 'baro'):
            from yaff.sampling.nvt import NHCThermostat
            if isinstance(self.thermostat, NHCThermostat):
                # in case the barostat is coupled with a NHC thermostat:
                # v_{xi,1} is needed to update v_g
                self.chainvel0 = self.thermostat.chain.vel[0]
            self.barostat.post(iterative, self.chainvel0)
        # update the correction on E_cons due to thermostat and barostat
        self.econs_correction = self.thermostat.econs_correction + self.barostat.econs_correction
        if isinstance(self.thermostat, NHCThermostat):
            # extra correction necessary if NHC thermostat is used
            kt = boltzmann*self.thermostat.temp
            self.econs_correction += self.barostat.dim**2*kt*self.thermostat.chain.pos[0]

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


class MartynaTobiasKleinBarostat(VerletHook):
    def __init__(self, ff, temp, press, start=0, timecon=1000*femtosecond):
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
        """
        self.temp = temp
        self.press = press
        self.timecon_press = timecon
        self.cell = ff.system.cell.rvecs.copy()
        self.dim = ff.system.cell.nvec
        # number of degrees of freedom is not yet known
        self.dof = 0
        self.mass_press = 0
        self.vel_press = np.zeros((3,3),float)

        # symmetrize the cell tensor
        self.cell_symmetrize(ff)

        VerletHook.__init__(self, start, 1)

    def set_ndof(self, ndof):
        # allocate degrees of freedom
        angfreq = 2*np.pi/self.timecon_press
        self.mass_press = (ndof+self.dim**2)*boltzmann*self.temp/angfreq**2
        self.vel_press = self.get_random_vel_press()

    def cell_symmetrize(self, ff):
        pos_old = ff.system.pos.copy()
        # frac_pos_old = np.zeros((len(pos_old),3), float)
        # frac_pos_new = np.zeros((len(pos_old),3), float)
        # for i in np.arange(0,len(pos_old)):
        #    frac_pos_old[i] = np.dot(np.linalg.inv(self.cell), pos_old[i])
        U, s, V = np.linalg.svd(self.cell)
        rot_mat = np.dot(V.T, U.T)
        self.cell = np.dot(rot_mat,self.cell)
        ff.update_rvecs(self.cell)
        pos_new = pos_old
        for i in np.arange(0,len(pos_old)):
            pos_new[i] = np.dot(rot_mat,pos_old[i])
        #    frac_pos_new[i] = np.dot(np.linalg.inv(self.cell), pos_new[i])
        ff.update_pos(pos_new)
        # print frac_pos_old - frac_pos_new

    def get_random_vel_press(self):
        # generates symmetric tensor of barostat velocities
        shape = 3, 3
        # generate random 3x3 tensor
        rand = np.random.normal(0, np.sqrt(self.mass_press*boltzmann*self.temp), shape)/self.mass_press
        vel_press = np.zeros(shape)
        # create initial symmetric pressure velocity tensor
        for i in np.arange(0,3):
            for j in np.arange(0,3):
                if i >= j:
                    vel_press[i,j] = rand[i,j]
                else:
                    vel_press[i,j] = rand[j,i]
        return vel_press

    def init(self, iterative):
        self.timestep_press = iterative.timestep

    def pre(self, iterative):
        pass

    def post(self, iterative):
        pass

    def propagate_press(self, chain_vel, ndof, ekin, vel, masses, volume, iterative):
        # iL vxi_1 h/8
        self.vel_press *= np.exp(-chain_vel*self.timestep_press/8)

        # necessary to calculate it here again, instead of during Verlet step?
        gpos = np.zeros(iterative.pos.shape, float)
        vtens = np.zeros((3,3),float)

        energy = iterative.ff.compute(gpos,vtens)
        #print str(iterative.vtens/kjmol)
        #ontbrekend = np.zeros((3,3),float)
        #for i in np.arange(0,len(iterative.pos)):
        #    ontbrekend += np.outer(-gpos[i],iterative.pos[i])
        ptens_vol = (np.dot(vel.T*masses, vel) - vtens)
        #print 'pxp = ' + str(np.dot(vel.T*masses,vel)/kjmol)
        #print 'vtens = ' + str(vtens/kjmol)
        ptens_vol = 0.5*(ptens_vol.T + ptens_vol)
        # print np.dot(vel.T*masses, vel)
        # print ptens_vol
        # print self.press*volume-2.0*ekin/ndof
        # print self.press*volume
        G = (ptens_vol+(2.0*ekin/ndof-self.press*volume)*np.eye(3))/self.mass_press
        #print 'PintV = ' + str(ptens_vol/kjmol)
        #print 'PextV - 2K/3N = ' + str((self.press*volume-2.0*ekin/ndof)/kjmol)
        #print 'ndof = ' + str(ndof)
        #print 'PextV = ' + str((self.press*volume)/kjmol)
        #print '2K/3N = ' + str((2.0*ekin/ndof)/kjmol)
        # print 'vol = ' + str(volume)
        # G = (ptens_vol-(self.press*volume)*np.eye(3))/self.mass_press
        # iL G_g h/4
        self.vel_press += G*self.timestep_press/4
        # iL vxi_1 h/8
        self.vel_press *= np.exp(-chain_vel*self.timestep_press/8)

    def propagate_vel(self, chain_vel, ndof, vel, masses):
        # diagonalize propagator matrix
        #print self.vel_press+(np.trace(self.vel_press)/ndof+chain_vel)*np.eye(3)
        Dg, Eg = np.linalg.eigh(self.vel_press+(np.trace(self.vel_press)/ndof+chain_vel)*np.eye(3))
        #Dg, Eg = np.linalg.eig(self.vel_press+(chain_vel)*np.eye(3))
        # define D_g and D'_g
        Daccg = np.exp(-Dg*self.timestep_press/2)
        Daccg = np.diagflat(Daccg)
        # iL (vg + Tr(vg)/ndof + vxi_1) h/2
        # and update kinetic energie
        ekin = 0
        for i in np.arange(0,len(vel)):
            vel[i] = np.dot(Eg, np.dot(Daccg, np.dot(Eg.T, vel[i])))
        ekin = 0.5*(vel**2*masses.reshape(-1,1)).sum()
        return vel, ekin

    def add_press_cont(self):
        # pressure contribution to g1: kinetic cell tensor energy
        # and extra degrees of freedom due to cell tensor
        return self.mass_press*np.trace(np.dot(self.vel_press.T,self.vel_press)) - self.dim**2*self.temp*boltzmann
        #return self.mass_press*np.trace(np.dot(self.vel_press.T,self.vel_press)) - self.dim*(self.dim+1)/2*self.temp*boltzmann

    def get_econs_correction(self, chain_pos, volume):
        kt = boltzmann*self.temp
        # add correction due to combination barostat and thermostat
        return self.dim**2*kt*chain_pos + 0.5*self.mass_press*np.trace(np.dot(self.vel_press.T,self.vel_press)) + self.press*volume
        #return self.dim*(self.dim+1)/2*kt*chain_vel + 0.5*self.mass_press*np.trace(np.dot(self.vel_press.T,self.vel_press)) + self.press*iterative.ff.system.cell.volume
