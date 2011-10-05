# YAFF is yet another force-field code
# Copyright (C) 2008 - 2011 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
# for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
# reserved unless otherwise stated.
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
# --


import numpy as np

from yaff import kjmol, angstrom, deg, angstrom
from yaff import ForceField, ForcePartValence

from yaff.pes.test.common import get_system_water32


def test_generator_water32_bondharm():
    system = get_system_water32()
    ff = ForceField.generate(system, 'input/parameters_water_bondharm.txt')
    assert len(ff.parts) == 1
    assert isinstance(ff.parts[0], ForcePartValence)
    part_valence = ff.parts[0]
    assert part_valence.dlist.ndelta == 64
    for i, j in system.topology.bonds:
        row = part_valence.dlist.lookup.get((i, j))
        assert row is not None
    assert (part_valence.iclist.ictab['kind'] == 0).all()
    assert part_valence.iclist.nic == 64
    assert (part_valence.vlist.vtab['kind'] == 0).all()
    assert abs(part_valence.vlist.vtab['par0'] - 4.0088096730e+03*(kjmol/angstrom**2)).max() < 1e-10
    assert abs(part_valence.vlist.vtab['par1'] - 1.0238240000e+00*angstrom).max() < 1e-10
    assert part_valence.vlist.nv == 64


def test_generator_water32_bondfues():
    system = get_system_water32()
    ff = ForceField.generate(system, 'input/parameters_water_bondfues.txt')
    assert len(ff.parts) == 1
    assert isinstance(ff.parts[0], ForcePartValence)
    part_valence = ff.parts[0]
    assert part_valence.dlist.ndelta == 64
    for i, j in system.topology.bonds:
        row = part_valence.dlist.lookup.get((i, j))
        assert row is not None
    assert (part_valence.iclist.ictab['kind'] == 0).all()
    assert part_valence.iclist.nic == 64
    assert (part_valence.vlist.vtab['kind'] == 2).all()
    assert abs(part_valence.vlist.vtab['par0'] - 4.0088096730e+03*(kjmol/angstrom**2)).max() < 1e-10
    assert abs(part_valence.vlist.vtab['par1'] - 1.0238240000e+00*angstrom).max() < 1e-10
    assert part_valence.vlist.nv == 64


def test_generator_water32_bendaharm():
    system = get_system_water32()
    ff = ForceField.generate(system, 'input/parameters_water_bendaharm.txt')
    assert len(ff.parts) == 1
    assert isinstance(ff.parts[0], ForcePartValence)
    part_valence = ff.parts[0]
    assert part_valence.dlist.ndelta == 64
    for i, j in system.topology.bonds:
        row = part_valence.dlist.lookup.get((i, j))
        assert row is not None
    assert (part_valence.iclist.ictab['kind'] == 2).all()
    assert part_valence.iclist.nic == 32
    assert (part_valence.vlist.vtab['kind'] == 0).all()
    assert abs(part_valence.vlist.vtab['par0'] - 3.0230353700e+02*kjmol).max() < 1e-10
    assert abs(part_valence.vlist.vtab['par1'] - 8.8401698835e+01*deg).max() < 1e-10
    assert part_valence.vlist.nv == 32


def test_generator_water32_bendcharm():
    system = get_system_water32()
    ff = ForceField.generate(system, 'input/parameters_water_bendcharm.txt')
    assert len(ff.parts) == 1
    assert isinstance(ff.parts[0], ForcePartValence)
    part_valence = ff.parts[0]
    assert part_valence.dlist.ndelta == 64
    for i, j in system.topology.bonds:
        row = part_valence.dlist.lookup.get((i, j))
        assert row is not None
    assert (part_valence.iclist.ictab['kind'] == 1).all()
    assert part_valence.iclist.nic == 32
    assert (part_valence.vlist.vtab['kind'] == 0).all()
    assert abs(part_valence.vlist.vtab['par0'] - 3.0230353700e+02*kjmol).max() < 1e-10
    assert abs(part_valence.vlist.vtab['par1'] - np.cos(8.8401698835e+01*deg)).max() < 1e-10
    assert part_valence.vlist.nv == 32


def test_generator_water32_fixq():
    system = get_system_water32()
    system.charges[:] = 0.0
    ff = ForceField.generate(system, 'input/parameters_water_fixq.txt', rcut=15.0*angstrom)
    assert len(ff.parts) == 4
    part_pair_ei = ff.part_pair_ei
    part_ewald_reci = ff.part_ewald_reci
    part_ewald_cor = ff.part_ewald_cor
    part_ewald_neut = ff.part_ewald_neut
    # check part settings
    assert part_pair_ei.pair_pot.alpha > 0
    assert part_pair_ei.pair_pot.alpha == part_ewald_reci.alpha
    assert part_pair_ei.pair_pot.alpha == part_ewald_cor.alpha
    assert part_pair_ei.pair_pot.alpha == part_ewald_neut.alpha
    # check charges
    for i in xrange(system.natom):
        if system.numbers[i] == 1:
            assert abs(system.charges[i] - 0.417) < 1e-5
        else:
            assert abs(system.charges[i] + 2*0.417) < 1e-5

    system = get_system_water32()
    system.charges[:] = 0.0
    ff2 = ForceField.generate(system, 'input/parameters_water_fixq.txt', rcut=15.0*angstrom)
    energy = ff.compute()
    energy2 = ff2.compute()
    print energy
    print energy2
    assert abs(energy - energy2) < 1e-3


def test_generator_water32_exprep1():
    system = get_system_water32()
    ff = ForceField.generate(system, 'input/parameters_water_exprep1.txt')
    assert len(ff.parts) == 1
    part_pair_exprep = ff.part_pair_exprep
    # check parameters
    assert part_pair_exprep.pair_pot.amp_mix == 0
    assert part_pair_exprep.pair_pot.amp_mix_coeff == 0.0
    assert part_pair_exprep.pair_pot.b_mix == 0
    assert part_pair_exprep.pair_pot.b_mix_coeff == 0.0
    assert part_pair_exprep.pair_pot.amps[0] == 4.2117588157e+02
    assert part_pair_exprep.pair_pot.amps[1] == 2.3514195495e+00
    assert abs(part_pair_exprep.pair_pot.bs[0] - 4.4661933834e+00/angstrom) < 1e-10
    assert abs(part_pair_exprep.pair_pot.bs[1] - 4.4107388814e+00/angstrom) < 1e-10


def test_generator_water32_exprep2():
    system = get_system_water32()
    ff = ForceField.generate(system, 'input/parameters_water_exprep2.txt')
    assert len(ff.parts) == 1
    part_pair_exprep = ff.part_pair_exprep
    # check parameters
    assert part_pair_exprep.pair_pot.amp_mix == 1
    assert part_pair_exprep.pair_pot.amp_mix_coeff == 2.385e-2
    assert part_pair_exprep.pair_pot.b_mix == 1
    assert part_pair_exprep.pair_pot.b_mix_coeff == 7.897e-3
    assert part_pair_exprep.pair_pot.amps[0] == 4.2117588157e+02
    assert part_pair_exprep.pair_pot.amps[1] == 2.3514195495e+00
    assert abs(part_pair_exprep.pair_pot.bs[0] - 4.4661933834e+00/angstrom) < 1e-10
    assert abs(part_pair_exprep.pair_pot.bs[1] - 4.4107388814e+00/angstrom) < 1e-10


def test_generator_water32():
    system = get_system_water32()
    system.charges[:] = 0.0
    ff = ForceField.generate(system, 'input/parameters_water.txt')
    # get all ff parts
    assert len(ff.parts) == 6
    part_valence = ff.part_valence
    part_pair_exprep = ff.part_pair_exprep
    part_pair_ei = ff.part_pair_ei
    part_ewald_reci = ff.part_ewald_reci
    part_ewald_cor = ff.part_ewald_cor
    part_ewald_neut = ff.part_ewald_neut
    # check exprep parameters
    assert part_pair_exprep.pair_pot.amp_mix == 1
    assert part_pair_exprep.pair_pot.amp_mix_coeff == 2.385e-2
    assert part_pair_exprep.pair_pot.b_mix == 1
    assert part_pair_exprep.pair_pot.b_mix_coeff == 7.897e-3
    assert part_pair_exprep.pair_pot.amps[0] == 4.2117588157e+02
    assert part_pair_exprep.pair_pot.amps[1] == 2.3514195495e+00
    assert abs(part_pair_exprep.pair_pot.bs[0] - 4.4661933834e+00/angstrom) < 1e-10
    assert abs(part_pair_exprep.pair_pot.bs[1] - 4.4107388814e+00/angstrom) < 1e-10
    # check charges
    for i in xrange(system.natom):
        if system.numbers[i] == 1:
            assert abs(system.charges[i] - 3.6841957737e+00) < 1e-5
        else:
            assert abs(system.charges[i] + 2*3.6841957737e+00) < 1e-5
    # check valence
    assert part_valence.dlist.ndelta == 64
    assert part_valence.iclist.nic == 96
    assert (part_valence.iclist.ictab['kind'][:64] == 0).all()
    assert (part_valence.iclist.ictab['kind'][64:96] == 1).all()
    assert part_valence.vlist.nv == 96
    assert abs(part_valence.vlist.vtab['par0'][:64] - 4.0088096730e+03*(kjmol/angstrom**2)).max() < 1e-10
    assert abs(part_valence.vlist.vtab['par1'][:64] - 1.0238240000e+00*angstrom).max() < 1e-10
    assert abs(part_valence.vlist.vtab['par0'][64:96] - 3.0230353700e+02*kjmol).max() < 1e-10
    assert abs(part_valence.vlist.vtab['par1'][64:96] - np.cos(8.8401698835e+01*deg)).max() < 1e-10
