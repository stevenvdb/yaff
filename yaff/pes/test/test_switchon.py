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


import numpy as np
from scipy.special import erf

from molmod.minimizer import check_delta

from yaff import *
from yaff.test.common import get_system_water32
from yaff.pes.test.common import check_gpos_part, check_vtens_part, check_hess_part
from test_pair_pot import check_pair_pot_water32

def get_part_water32_lj_switcherf():
    # Initialize system, nlist and scaling
    system = get_system_water32()
    nlist = NeighborList(system)
    scalings = Scalings(system)
    # Set atomic radii
    system.radii = np.zeros(system.natom)
    system.radii[system.numbers==1] = 0.7*angstrom
    system.radii[system.numbers==8] = 1.4*angstrom
    # Initialize parameters
    rminhalf_table = {1: 0.2245*angstrom, 8: 1.7682*angstrom}
    epsilon_table = {1: -0.0460*kcalmol, 8: -0.1521*kcalmol}
    sigmas = np.zeros(96, float)
    epsilons = np.zeros(96, float)
    for i in xrange(system.natom):
        sigmas[i] = rminhalf_table[system.numbers[i]]*(2.0)**(5.0/6.0)
        epsilons[i] = epsilon_table[system.numbers[i]]
    # Create the pair_pot and part_pair
    rcut = 9*angstrom
    sw = SwitchErrorFunction(system.radii)
    pair_pot = PairPotLJ(sigmas, epsilons, rcut, Hammer(1.0), sw=sw)
    assert abs(pair_pot.sigmas - sigmas).max() == 0.0
    assert abs(pair_pot.epsilons - epsilons).max() == 0.0
    part_pair = ForcePartPair(system, nlist, scalings, pair_pot)
    # Create a pair function:
    def pair_fn(i, j, d, delta):
        sigma = 0.5*(sigmas[i]+sigmas[j])
        epsilon = np.sqrt(epsilons[i]*epsilons[j])
        x = (sigma/d)**6
        radius = np.sqrt(system.radii[i]**2 + system.radii[j]**2)
        return 4*epsilon*(x*(x-1))*np.exp(1.0/(d-rcut))*erf(d/radius)
    return system, nlist, scalings, part_pair, pair_fn

def test_switchon_water32():
    system, nlist, scalings, part_pair, pair_fn = get_part_water32_lj_switcherf()
    check_pair_pot_water32(system, nlist, scalings, part_pair, pair_fn, 1e-15)
    check_gpos_part(system, part_pair, nlist)
    check_vtens_part(system, part_pair, nlist)
    check_hess_part(system, part_pair, nlist)
