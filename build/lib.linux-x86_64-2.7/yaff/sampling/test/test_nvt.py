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


from yaff import *
from yaff.sampling.test.common import get_ff_water32


def test_basic_nhc_nvt():
    thermostat = NHCThermostat(300)
    nvt = VerletIntegrator(get_ff_water32(), 1.0*femtosecond, hooks=thermostat)
    nvt.run(5)
    assert nvt.counter == 5

def test_basic_langevin_nvt():
    thermostat = LangevinThermostat(300)
    nvt = VerletIntegrator(get_ff_water32(), 1.0*femtosecond, hooks=thermostat)
    nvt.run(5)
    assert nvt.counter == 5

def test_at():
    nve = VerletIntegrator(get_ff_water32(), 1.0*femtosecond, hooks=AndersenThermostat(300))
    nve.run(5)
    assert nve.counter == 5


def test_at_select():
    nve = VerletIntegrator(get_ff_water32(), 1.0*femtosecond, hooks=AndersenThermostat(300, select=[1,2,5]))
    nve.run(5)
    assert nve.counter == 5


def test_at_annealing():
    nve = VerletIntegrator(get_ff_water32(), 1.0*femtosecond, hooks=AndersenThermostat(300, annealing=0.9))
    nve.run(5)
    assert nve.counter == 5
