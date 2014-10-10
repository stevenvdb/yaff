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
'''Thermostats'''


import numpy as np

from molmod import boltzmann, femtosecond

from yaff.log import log
from yaff.sampling.iterative import Iterative, StateItem
from yaff.sampling.utils import get_random_vel, clean_momenta, \
    get_ndof_internal_md
from yaff.sampling.verlet import VerletHook


__all__ = [
    'AndersenThermostat', 'BerendsenThermostat', 'LangevinThermostat',
    'CSVRThermostat', 'NHCThermostat', 'NHCAttributeStateItem',
]


class AndersenThermostat(VerletHook):
    def __init__(self, temp, start=0, step=1, select=None, annealing=1.0):
        """
           This is an implementation of the Andersen thermostat. The method
           is described in:

                Andersen, H. C. J. Chem. Phys. 1980, 72, 2384-2393.

           **Arguments:**

           temp
                The average temperature of the NVT ensemble

           **Optional arguments:**

           start
                The first iteration at which this hook is called

           step
                The number of iterations between two subsequent calls to this
                hook.

           select
                An array of atom indexes to indicate which atoms controlled by
                the thermostat.

           annealing
                After every call to this hook, the temperature is multiplied
                with this annealing factor. This effectively cools down the
                system.
        """
        self.temp = temp
        self.select = select
        self.annealing = annealing
        VerletHook.__init__(self, start, step)

    def init(self, iterative):
        # It is mandatory to zero the external momenta.
        clean_momenta(iterative.pos, iterative.vel, iterative.masses, iterative.ff.system.cell)

    def pre(self, iterative, G1_add = None):
        # Andersen thermostat step before usual Verlet hook, since it largely affects the velocities
        # Needed to correct the conserved quantity
        ekin_before = iterative._compute_ekin()
        # Change the (selected) velocities
        if self.select is None:
            iterative.vel[:] = get_random_vel(self.temp, False, iterative.masses)
        else:
            iterative.vel[self.select] = get_random_vel(self.temp, False, iterative.masses, self.select)
        # Zero any external momenta after choosing new velocities
        clean_momenta(iterative.pos, iterative.vel, iterative.masses, iterative.ff.system.cell)
        # Update the kinetic energy and the reference for the conserved quantity
        ekin_after = iterative._compute_ekin()
        self.econs_correction += ekin_before - ekin_after
        # Optional annealing
        self.temp *= self.annealing

    def post(self, iterative, G1_add = None):
        pass

class BerendsenThermostat(VerletHook):
    def __init__(self, temp, start=0, timecon=100*femtosecond):
        """
           This is an implementation of the Berendsen thermostat. The algorithm
           is described in:

                Berendsen, H. J. C.; Postma, J. P. M.; van Gunsteren, W. F.;
                Dinola, A.; Haak, J. R. J. Chem. Phys. 1984, 81, 3684-3690

           **Arguments:**

           temp
                The temperature of thermostat.

           **Optional arguments:**

           start
                The step at which the thermostat becomes active.

           timecon
                The time constant of the Berendsen thermostat.
        """
        self.temp = temp
        self.timecon = timecon
        VerletHook.__init__(self, start, 1)

    def init(self, iterative):
        # It is mandatory to zero the external momenta.
        clean_momenta(iterative.pos, iterative.vel, iterative.masses, iterative.ff.system.cell)
        if iterative.ndof is None:
            iterative.ndof = get_ndof_internal_md(iterative.pos.shape[0], iterative.ff.system.cell.nvec)

    def pre(self, iterative, G1_add = None):
        ekin = iterative.ekin
        temp_inst = 2.0*iterative.ekin/(boltzmann*iterative.ndof)
        c = np.sqrt(1+iterative.timestep/self.timecon*(self.temp/temp_inst-1))
        iterative.vel[:] = c*iterative.vel
        iterative.ekin = iterative._compute_ekin()
        self.econs_correction += (1-c**2)*ekin

    def post(self, iterative, G1_add = None):
        pass

class LangevinThermostat(VerletHook):
    def __init__(self, temp, start=0, timecon=100*femtosecond):
        """
           This is an implementation of the Langevin thermostat. The algorithm
           is described in:

                Bussi, G.; Parrinello, M. Phys. Rev. E 2007, 75, 056707

           **Arguments:**

           temp
                The temperature of thermostat.

           **Optional arguments:**

           start
                The step at which the thermostat becomes active.

           timecon
                The time constant of the Nose-Hoover thermostat.
        """
        self.temp = temp
        self.timecon = timecon
        VerletHook.__init__(self, start, 1)

    def init(self, iterative):
        # It is mandatory to zero the external momenta.
        clean_momenta(iterative.pos, iterative.vel, iterative.masses, iterative.ff.system.cell)
        # Necessary to calculate conserved quantity
        #self.forces = np.zeros(iterative.pos.shape, float)
        #self.posoud = iterative.pos.copy()

    def pre(self, iterative, G1_add = None):
        # Necessary to calculate conserved quantity
        #self.epot_old = iterative.epot
        #self.forces[:] = iterative.acc*iterative.masses.reshape(-1,1)
        ekin0 = iterative.ekin
        # Actual update
        self.thermo(iterative)
        ekin1 = iterative.ekin
        iterative.ekin += ekin0-ekin1

    def post(self, iterative, G1_add = None):
        ekin0 = iterative.ekin
        # Actual update
        self.thermo(iterative)
        ekin1 = iterative.ekin
        iterative.ekin += ekin0-ekin1
        # Calculation of the conserved quantity
        #forces_new = iterative.acc*iterative.masses.reshape(-1,1)
        #self.econs_correction = -iterative.ekin + iterative.econs + 0.5*np.multiply(iterative.pos.copy()-self.posoud, self.forces + forces_new).sum() - self.epot_old + iterative.timestep**2/8.0*np.divide(np.sum(forces_new**2-self.forces**2, axis=1), iterative.masses).sum()
        #self.posoud = iterative.pos.copy()

    def thermo(self, iterative):
        c1 = np.exp(-iterative.timestep/self.timecon/2)
        c2 = np.sqrt((1.0-c1**2)*self.temp*boltzmann/iterative.masses).reshape(-1,1)
        iterative.vel[:] = c1*iterative.vel + c2*np.random.normal(0, 1, iterative.vel.shape)
        iterative.ekin = iterative._compute_ekin()

class CSVRThermostat(VerletHook):
    def __init__(self, temp, start=0, timecon=100*femtosecond):
        """
            This is an implementation of the CSVR thermostat. The equations are
            derived in:

                Bussi, G.; Donadio, D.; Parrinello, M. J. Chem. Phys. 2007,
                126, 014101

            The implementation (used here) is derived in

                Bussi, G.; Parrinello, M. Comput. Phys. Commun. 2008, 179, 26-29

           **Arguments:**

           temp
                The temperature of thermostat.

           **Optional arguments:**

           start
                The step at which the thermostat becomes active.

           timecon
                The time constant of the CSVR thermostat.
        """
        self.temp = temp
        self.timecon = timecon
        VerletHook.__init__(self, start, 1)

    def init(self, iterative):
        # It is mandatory to zero the external momenta.
        clean_momenta(iterative.pos, iterative.vel, iterative.masses, iterative.ff.system.cell)
        if iterative.ndof is None:
            iterative.ndof = get_ndof_internal_md(iterative.pos.shape[0], iterative.ff.system.cell.nvec)
        self.kin = 0.5*iterative.ndof*boltzmann*self.temp

    def pre(self, iterative, G1_add = None):
        c = np.exp(-iterative.timestep/self.timecon)
        R = np.random.normal(0, 1)
        S = (np.random.normal(0, 1, iterative.ndof-1)**2).sum()
        iterative.ekin = iterative._compute_ekin()
        fact = (1-c)*self.kin/iterative.ndof/iterative.ekin
        alpha = np.sign(R+np.sqrt(c/fact))*np.sqrt(c + (S+R**2)*fact + 2*R*np.sqrt(c*fact))
        iterative.vel[:] = alpha*iterative.vel
        iterative.ekin_new = alpha**2*iterative.ekin
        self.econs_correction += (1-alpha**2)*iterative.ekin
        iterative.ekin = iterative.ekin_new

    def post(self, iterative, G1_add = None):
        pass
        # How implementing the Wiener Noise?
        # self.econs_correction += (self.kin-iterative.ekin)*iterative.timestep/self.timecon


class NHChain(object):
    def __init__(self, length, timestep, temp, ndof, timecon=100*femtosecond):
        # parameters
        self.length = length
        self.timestep = timestep
        self.temp = temp
        self.timecon = timecon
        if ndof>0:#avoid setting self.masses with zero gaussian-width in set_ndof if ndof=0
            self.set_ndof(ndof)

        # allocate degrees of freedom
        self.pos = np.zeros(length)
        self.vel = np.zeros(length)

    def set_ndof(self, ndof):
        # set the masses according to the time constant
        self.ndof = ndof
        angfreq = 2*np.pi/self.timecon
        self.masses = np.ones(self.length)*(boltzmann*self.temp/angfreq**2)
        self.masses[0] *= ndof
        self.vel = self.get_random_vel_therm()

    def get_random_vel_therm(self):
        # generate random velocities for the thermostat velocities using a Gaussian distribution
        shape = self.length
        return np.random.normal(0, np.sqrt(self.masses*boltzmann*self.temp), shape)/self.masses

    def __call__(self, ekin, vel, G1_add):
        def do_bead(k, ekin):
            # Compute g
            if k == 0:
                # coupling with atoms (and barostat)
                # L = ndof (+d(d+1)/2 (aniso NPT) / +1 (iso NPT)) because of equidistant time steps.
                g = 2*ekin - self.ndof*self.temp*boltzmann
                if G1_add is not None:
                    # add pressure contribution to g1
                    g += G1_add
            else:
                # coupling between beads
                g = self.masses[k-1]*self.vel[k-1]**2 - self.temp*boltzmann
            g /= self.masses[k]

            # Lioville operators on relevant part of the chain
            if k == self.length-1:
                # iL G_k h/4
                self.vel[k] += g*self.timestep/4
            else:
                # iL vxi_{k-1} h/8
                self.vel[k] *= np.exp(-self.vel[k+1]*self.timestep/8)
                # iL G_k h/4
                self.vel[k] += g*self.timestep/4
                # iL vxi_{k-1} h/8
                self.vel[k] *= np.exp(-self.vel[k+1]*self.timestep/8)

        # Loop over chain in reverse order
        for k in xrange(self.length-1, -1, -1):
            do_bead(k, ekin)

        # iL xi (all) h/2
        self.pos += self.vel*self.timestep/2
        # iL Cv (all) h/2
        factor = np.exp(-self.vel[0]*self.timestep/2)
        vel *= factor
        ekin *= factor**2

        # Loop over chain in forward order
        for k in xrange(0, self.length):
            do_bead(k, ekin)
        return vel, ekin

    def get_econs_correction(self):
        kt = boltzmann*self.temp
        # correction due to the thermostat
        return 0.5*(self.vel**2*self.masses).sum() + kt*(self.ndof*self.pos[0] + self.pos[1:].sum())


class NHCThermostat(VerletHook):
    def __init__(self, temp, start=0, timecon=100*femtosecond, chainlength=3):
        """
           This hook implements the Nose-Hoover-Chain thermostat. The equations
           are derived in:

                Martyna, G. J.; Klein, M. L.; Tuckerman, M. J. Chem. Phys. 1992,
                97, 2635-2643.

           The implementation (used here) of a symplectic integrator of the
           Nose-Hoover-Chain thermostat is discussed in:

                Martyna, G. J.;  Tuckerman, M. E.;  Tobias, D. J.;  Klein,
                M. L. Mol. Phys. 1996, 87, 1117-1157.

           **Arguments:**

           temp
                The temperature of thermostat.

           **Optional arguments:**

           start
                The step at which the thermostat becomes active.

           timecon
                The time constant of the Nose-Hoover thermostat.

           chainlength
                The number of beads in the Nose-Hoover chain.

        """
        self.temp = temp
        # At this point, the timestep and the number of degrees of freedom are
        # not known yet.
        self.chain = NHChain(chainlength, 0.0, temp, 0, timecon)
        VerletHook.__init__(self, start, 1)

    def init(self, iterative):
        # It is mandatory to zero the external momenta.
        clean_momenta(iterative.pos, iterative.vel, iterative.masses, iterative.ff.system.cell)
        # If needed, determine the number of _internal_ degrees of freedom.
        if iterative.ndof is None:
            iterative.ndof = get_ndof_internal_md(iterative.pos.shape[0], iterative.ff.system.cell.nvec)
        # Configure the chain.
        self.chain.timestep = iterative.timestep
        self.chain.set_ndof(iterative.ndof)

    def pre(self, iterative, G1_add = None):
        vel_new, iterative.ekin = self.chain(iterative.ekin, iterative.vel, G1_add)
        iterative.vel[:] = vel_new

    def post(self, iterative, G1_add = None):
        vel_new, iterative.ekin = self.chain(iterative.ekin, iterative.vel, G1_add)
        iterative.vel[:] = vel_new
        self.econs_correction = self.chain.get_econs_correction()

class NHCAttributeStateItem(StateItem):
    def __init__(self, attr):
        StateItem.__init__(self, 'chain_'+attr)
        self.attr = attr

    def get_value(self, iterative):
        chain = None
        for hook in iterative.hooks:
            if isinstance(hook, NHCThermostat):
                chain = hook.chain
                break
            elif isinstance(hook, TBCombination):
                if isinstance(hook.thermo, NHCThermostat):
                    chain = hook.thermo.chain
                break
        if chain is None:
            raise TypeError('Iterative does not contain a NHCThermostat hook.')
        return getattr(chain, self.attr)
