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

from molmod import kcalmol, angstrom, rad, deg, femtosecond, boltzmann

from yaff.test.common import get_system_water_polarizable
from yaff.sampling.polarization import Polarization_FSPM_ISD, construct_ei_tensor_monomono, \
        construct_ei_tensor_dipdip, construct_ei_tensor_monodip
from yaff.pes import *
from yaff.sampling.iterative import Iterative

from yaff.pes.test.test_pair_pot import get_part_water_eidip
from yaff.sampling.polarization import *

def get_part_water_eidip(scalings = [1.0,1.0,1.0],rcut=14.0*angstrom,switch_width=0.0*angstrom, alpha=0.0):
    # Initialize system, nlist and scaling
    system = get_system_water_polarizable()
    nlist = NeighborList(system)
    scalings = Scalings(system, scalings[0], scalings[1], scalings[2])
    # Set dipoles
    slater1p_N = np.zeros( (system.natom,3) )
    slater1p_Z = np.zeros( (system.natom,3) )
    polarizabilities = np.array([6.6,1.6,1.6])
    slater1p_widths = np.zeros( (system.natom,3) ) # This should be related to polarizability
    for i in xrange(system.natom):
        slater1p_widths[i,:] = np.power( 7.0/32.0/24.0*polarizabilities[i] , 1.0/3.0 )
    # Create the pair_pot and part_pair
    pair_pot_ei = PairPotEIDip(system.slater1s_Z+system.slater1s_N, slater1p_N+slater1p_Z, alpha, rcut, tr = Switch3(switch_width))
    part_pair_ei = ForcePartPair(system, nlist, scalings, pair_pot_ei)
    # Pair potential that contains the Slater corrections
    pair_pot_slatercorr = PairPotEiSlater1sp1spCorr(system.slater1s_widths,system.slater1s_N,system.slater1s_Z,
                                                    slater1p_widths, slater1p_N, slater1p_Z, rcut, tr = Switch3(switch_width))
    part_pair_slatercorr = ForcePartPair(system, nlist, scalings, pair_pot_slatercorr)
    nlist.update()
    return system, nlist, scalings, part_pair_ei, part_pair_slatercorr


def test_polarization_FSPM_ISD():
    '''
    Check that induced dipoles are set from our polarization model.
    '''
    system, nlist, scalings, part_pair_ei, part_pair_slatercorr = get_part_water_eidip()
    part_polarization = ForcePartPolarizationEnergy(system, None, None, None)
    ff = ForceField(system, [part_pair_ei,part_pair_slatercorr,part_polarization], nlist)
    iterative = Iterative(ff,hooks=Polarization_FSPM_ISD())
    # Check that at least one of the induced dipoles now has a non-zero value
    # We don't even check whether these values make any sense.
    assert np.any( part_pair_ei.pair_pot.dipoles != 0.0 )
    assert np.any( part_pair_slatercorr.pair_pot.slater1p_N != 0.0 )
    

def test_ei_tensor_monomono():
    system, nlist, scalings, part_pair_ei, part_pair_slatercorr = get_part_water_eidip()
    part_polarization = ForcePartPolarizationEnergy(system, None, None, None)
    ff = ForceField(system, [part_pair_ei,part_pair_slatercorr,part_polarization], nlist)
    iterative = Iterative(ff,hooks=Polarization_FSPM_ISD())
    iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1p_N[:] *= 0.0
    iterative.ff.part_pair_eidip.pair_pot.dipoles[:] *= 0.0
    # Compute the tensor
    Tff = construct_ei_tensor_monomono(iterative)
    # Setup vector with fixed multipoles
    Qf = np.zeros((2*system.natom,))
    Qf[::2] = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_Z.copy()
    Qf[1::2] = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_N.copy()
    # Compute energy using Yaff implementation in pes module
    Eff_check = iterative.ff.part_pair_eislater1sp1spcorr.compute() + iterative.ff.part_pair_eidip.compute()
    # Compute same energy using the tensor
    Eff_tensor = 0.5*np.dot(np.dot(Qf.transpose(),Tff),Qf)
    assert np.abs(Eff_check - Eff_tensor) < 1e-8


def test_ei_tensor_dipdip():
    system, nlist, scalings, part_pair_ei, part_pair_slatercorr = get_part_water_eidip()
    part_polarization = ForcePartPolarizationEnergy(system, None, None, None)
    ff = ForceField(system, [part_pair_ei,part_pair_slatercorr,part_polarization], nlist)
    iterative = Iterative(ff,hooks=Polarization_FSPM_ISD())
    iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_N[:] *= 0.0
    iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_Z[:] *= 0.0
    iterative.ff.part_pair_eidip.pair_pot.charges[:] *= 0.0
    # Compute the tensor
    Tff = construct_ei_tensor_dipdip(iterative)
    # Setup vector with fixed multipoles
    Qf = np.zeros((3*system.natom,))
    Qf[0::3] = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1p_N[:,0]
    Qf[1::3] = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1p_N[:,1]
    Qf[2::3] = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1p_N[:,2]
    # Compute energy using Yaff implementation in pes module
    Eff_check = iterative.ff.part_pair_eislater1sp1spcorr.compute() + iterative.ff.part_pair_eidip.compute()
    # Compute same energy using the tensor
    Eff_tensor = 0.5*np.dot(np.dot(Qf.transpose(),Tff),Qf)
    assert np.abs(Eff_check - Eff_tensor) < 1e-8


def test_ei_tensor_monodip():
    system, nlist, scalings, part_pair_ei, part_pair_slatercorr = get_part_water_eidip()
    part_polarization = ForcePartPolarizationEnergy(system, None, None, None)
    ff = ForceField(system, [part_pair_ei,part_pair_slatercorr,part_polarization], nlist)
    iterative = Iterative(ff,hooks=Polarization_FSPM_ISD())
    # Setup vector with monopoles
    Qf = np.zeros((2*system.natom,))
    Qf[::2] = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_Z.copy()
    Qf[1::2] = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_N.copy()
    # Setup vector with dipoles
    Qv = np.zeros((3*system.natom,))
    Qv[0::3] = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1p_N[:,0].copy()
    Qv[1::3] = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1p_N[:,1].copy()
    Qv[2::3] = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1p_N[:,2].copy()
    # Total interaction
    Etot = iterative.ff.part_pair_eislater1sp1spcorr.compute() + iterative.ff.part_pair_eidip.compute()
    # Dipole Dipole interaction
    iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_N[:] *= 0.0
    iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_Z[:] *= 0.0
    iterative.ff.part_pair_eidip.pair_pot.charges[:] *= 0.0
    Edd = iterative.ff.part_pair_eislater1sp1spcorr.compute() + iterative.ff.part_pair_eidip.compute()
    # Monopole Monopole interaction
    iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1p_N[:] *= 0.0
    iterative.ff.part_pair_eidip.pair_pot.dipoles[:] *= 0.0
    iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_N[:] = Qf[1::2]
    iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_Z[:] = Qf[0::2]
    iterative.ff.part_pair_eidip.pair_pot.charges[:] = Qf[0::2] + Qf[1::2]
    Eqq = iterative.ff.part_pair_eislater1sp1spcorr.compute() + iterative.ff.part_pair_eidip.compute()
    # Monopole Dipole interaction
    Eqd = Etot - Edd - Eqq
    # Calculate using tensor
    Tfv = construct_ei_tensor_monodip(iterative)
    E_tensor = 1.0*np.dot(np.dot(Qf.transpose(),Tfv),Qv)
    assert np.abs( Eqd - E_tensor ) < 1e-8


def test_polarization_DipolSCPicard():
    system, nlist, scalings, part_pair, pair_pot, pair_fn = get_part_water_eidip(scalings=[1.0,1.0,1.0])
    poltens_i = np.diag([1.0]*3*system.natom) #This is not used for this test
    dipoles = DipolSCPicard(system.pos, system.charges, poltens_i, system.natom, init=None, conv_crit=1e-10, tensors=None, system=system)
    G_0, G_1, G_2, D, chi = get_ei_tensors( system.pos, poltens_i, system.natom, system)
    lhs = np.dot( poltens_i + G_2,  np.reshape( dipoles , (-1) ) ) + chi
    rhs = np.dot(-G_1, system.charges)
    assert np.all( np.abs(lhs-rhs) < 1.0e-10 )

def test_tmp():
    #This is not really a test yet, just check if everything runs
    system, nlist, scalings, part_pair, pair_pot, pair_fn = get_part_water_eidip(scalings=[0.0,1.0,1.0])
    system.dipoles = np.zeros( (system.natom,3) )
    ff = ForceField(system, [part_pair], nlist)
    opt = CGOptimizer(CartesianDOF(ff), hooks=RelaxDipoles())
    opt.run(2)

def test_polarization_get_ei_tensors():
    """Check if the tensors from polarization module give correct energy"""
    #Don't scale interactions, this is not implemented in determining the tensors
    system, nlist, scalings, part_pair, pair_pot, pair_fn = get_part_water_eidip(scalings=[1.0,1.0,1.0])
    poltens_i = np.diag([1.0]*3*system.natom) #This is not used for this test
    #Get tensors from polarization module
    G_0, G_1, G_2, D, chi = get_ei_tensors( system.pos, poltens_i, system.natom, system)
    #Reshape the dipole matrix to simplify matrix expressions
    dipoles = np.reshape( pair_pot.dipoles , (-1,) )
    #Compute energy using these tensors
    #Charge-charge interaction
    energy_tensor = 0.5*np.dot(np.transpose(system.charges), np.dot(G_0,system.charges) )
    #Charge-dipole interaction
    energy_tensor += np.dot( np.transpose( dipoles), np.dot( G_1, system.charges) )
    #Dipole-dipole interaction
    energy_tensor += 0.5*np.dot( np.transpose( dipoles), np.dot( G_2, dipoles) )
    #Dipole creation energy
    energy_tensor += 0.5*np.dot( np.transpose( dipoles), np.dot( np.linalg.inv(D), dipoles) )
    nlist.update() # update the neighborlists, once the rcuts are known.
    # Compute the energy using yaff.
    energy_yaff = part_pair.compute()
    assert np.abs(energy_yaff-energy_tensor) < 1.0e-10
>>>>>>> Tests for polarization
