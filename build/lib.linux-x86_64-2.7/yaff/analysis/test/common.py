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


import tempfile, h5py as h5

from yaff import *
from yaff.sampling.test.common import get_ff_water32


def get_nve_water32():
    # Make a temporary directory
    dn_tmp = tempfile.mkdtemp(suffix='yaff', prefix='nve_water_32')
    # Setup a test FF
    ff = get_ff_water32()
    # Run a test simulation
    f = h5.File('%s/output.h5' % dn_tmp)
    hdf5 = HDF5Writer(f)
    nve = VerletIntegrator(ff, 1.0*femtosecond, hooks=hdf5)
    nve.run(5)
    assert nve.counter == 5
    return dn_tmp, nve, hdf5.f


def get_nvt_water32():
    # Make a temporary directory
    dn_tmp = tempfile.mkdtemp(suffix='yaff', prefix='nvt_water_32')
    # Setup a test FF
    ff = get_ff_water32()
    # Run a test simulation
    f = h5.File('%s/output.h5' % dn_tmp)
    hdf5 = HDF5Writer(f)
    thermostat = LangevinThermostat(temp=300)
    nvt = VerletIntegrator(ff, 1.0*femtosecond, hooks=[hdf5, thermostat])
    nvt.run(5)
    assert nvt.counter == 5
    return dn_tmp, nvt, hdf5.f


def get_opt_water32():
    # Make a temporary directory
    dn_tmp = tempfile.mkdtemp(suffix='yaff', prefix='opt_water_32')
    # Setup a test FF
    ff = get_ff_water32()
    # Run a test simulation
    f = h5.File('%s/output.h5' % dn_tmp)
    hdf5 = HDF5Writer(f)
    opt = CGOptimizer(FullCellDOF(ff), hooks=hdf5)
    opt.run(5)
    assert opt.counter == 5
    return dn_tmp, opt, hdf5.f
