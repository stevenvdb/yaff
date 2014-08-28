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
'''Polarizable forcefields
This is serious Work-In-Progress.
Very little attention is paid to efficient code, I am still in a phase of
testing out models...
'''
 
import numpy as np

from yaff.pes.ext import PairPotEIDip, PairPotEiSlater1sp1spCorr
from yaff.pes.ff import ForcePartPair
from yaff.sampling import Hook

__all__ = ['Polarization_FSPM_ISD']


def construct_ei_tensor_monomono(iterative):
    '''
    Construct the electrostatic interaction tensor with as basis functions:
        * a point monopole
        * a Slater monopole
    '''
    system = iterative.ff.system
    nbf = 2
    Tff = np.zeros((system.natom*nbf,system.natom*nbf)) # Tensor between fixed multipoles
    # Make some fake multipoles that help in calculating matrix elements
    # of the interaction tensor
    q1 = np.zeros((system.natom,))
    q2 = np.zeros((system.natom,))
    q3 = np.zeros((system.natom,))
    d1 = np.zeros((system.natom,3))
    d2 = np.zeros((system.natom,3))
    d3 = np.zeros((system.natom,3))
    # Construct potentials contributing to the tensor
    # Recycle settings from the force field
    alpha    = iterative.ff.part_pair_eidip.pair_pot.alpha
    rcut     = iterative.ff.part_pair_eidip.pair_pot.rcut
    tr       = iterative.ff.part_pair_eidip.pair_pot.get_truncation()
    nlist    = iterative.ff.part_pair_eidip.nlist
    scalings = iterative.ff.part_pair_eidip.scalings
    pair_pot_point  = PairPotEIDip(q1, d1, alpha, rcut, tr=tr )
    part_pair_point = ForcePartPair(system,  nlist, scalings, pair_pot_point)

    slater1s_widths = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_widths.copy()
    slater1p_widths = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1p_widths.copy()
    rcut            = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.rcut
    tr              = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.get_truncation()
    nlist           = iterative.ff.part_pair_eislater1sp1spcorr.nlist
    scalings        = iterative.ff.part_pair_eislater1sp1spcorr.scalings
    pair_pot_slatercorr  = PairPotEiSlater1sp1spCorr(slater1s_widths, q2, q3, slater1p_widths, d2, d3, rcut, tr=tr)
    part_pair_slatercorr = ForcePartPair(system, nlist, scalings, pair_pot_slatercorr)

    # Outer atom loop
    for a0 in xrange(system.natom):
        # Inner atom loop; because of symmetry we only consider a1>a0.
        # The self-interactions for the case a1=a0 are handled later on.
        for a1 in xrange(a0+1,system.natom):
            #
            # Construct fixed-fixed tensor for this atom pair
            #
            tff = np.zeros((nbf,nbf))
            # Point charge - Point charge
            pair_pot_point.charges[a0] = 1.0
            pair_pot_slatercorr.slater1s_Z[a0] = 1.0
            pair_pot_point.charges[a1] = 1.0
            pair_pot_slatercorr.slater1s_Z[a1] = 1.0
            tff[0,0] = part_pair_point.compute() + part_pair_slatercorr.compute()
            pair_pot_point.charges[:]  *= 0.0
            pair_pot_slatercorr.slater1s_Z[:] *= 0.0
            # Point charge - Slater charge
            pair_pot_point.charges[a0] = 1.0
            pair_pot_slatercorr.slater1s_Z[a0] = 1.0
            pair_pot_point.charges[a1] = 1.0
            pair_pot_slatercorr.slater1s_N[a1] = 1.0
            tff[0,1] = part_pair_point.compute() + part_pair_slatercorr.compute()
            pair_pot_point.charges[:]  *=0.0
            pair_pot_slatercorr.slater1s_Z[:] *= 0.0
            pair_pot_slatercorr.slater1s_N[:] *= 0.0
            # Slater charge - point charge
            pair_pot_point.charges[a0] = 1.0
            pair_pot_slatercorr.slater1s_N[a0] = 1.0
            pair_pot_point.charges[a1] = 1.0
            pair_pot_slatercorr.slater1s_Z[a1] = 1.0
            tff[1,0] = part_pair_point.compute() + part_pair_slatercorr.compute()
            pair_pot_point.charges[:]  *=0.0
            pair_pot_slatercorr.slater1s_Z[:] *= 0.0
            pair_pot_slatercorr.slater1s_N[:] *= 0.0
            # Slater charge - Slater charge
            pair_pot_point.charges[a0] = 1.0
            pair_pot_slatercorr.slater1s_N[a0] = 1.0
            pair_pot_point.charges[a1] = 1.0
            pair_pot_slatercorr.slater1s_N[a1] = 1.0
            tff[1,1] = part_pair_point.compute() + part_pair_slatercorr.compute()
            pair_pot_point.charges[:]  *=0.0
            pair_pot_slatercorr.slater1s_N[:] *= 0.0
            # Fill big tensor
            Tff[nbf*a0:nbf*(a0+1),nbf*a1:nbf*(a1+1)] = tff.copy()
            Tff[nbf*a1:nbf*(a1+1),nbf*a0:nbf*(a0+1)] = tff.transpose().copy()
    return Tff


def construct_ei_tensor_dipdip(iterative):
    '''
    Construct the electrostatic interaction tensor with as basis functions:
        * Slater dipoles with x, y and z orientation
    '''
    system = iterative.ff.system
    nbv = 3
    Tvv = np.zeros((system.natom*nbv,system.natom*nbv)) # Tensor between fixed multipoles
    # Make some fake multipoles that help in calculating matrix elements
    # of the interaction tensor
    q1 = np.zeros((system.natom,))
    q2 = np.zeros((system.natom,))
    q3 = np.zeros((system.natom,))
    d1 = np.zeros((system.natom,3))
    d2 = np.zeros((system.natom,3))
    d3 = np.zeros((system.natom,3))
    # Construct potentials contributing to the tensor
    # Recycle settings from the force field
    alpha    = iterative.ff.part_pair_eidip.pair_pot.alpha
    rcut     = iterative.ff.part_pair_eidip.pair_pot.rcut
    tr       = iterative.ff.part_pair_eidip.pair_pot.get_truncation()
    nlist    = iterative.ff.part_pair_eidip.nlist
    scalings = iterative.ff.part_pair_eidip.scalings
    pair_pot_point  = PairPotEIDip(q1, d1, alpha, rcut, tr=tr )
    part_pair_point = ForcePartPair(system,  nlist, scalings, pair_pot_point)

    slater1s_widths = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_widths.copy()
    slater1p_widths = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1p_widths.copy()
    rcut            = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.rcut
    tr              = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.get_truncation()
    nlist           = iterative.ff.part_pair_eislater1sp1spcorr.nlist
    scalings        = iterative.ff.part_pair_eislater1sp1spcorr.scalings
    pair_pot_slatercorr  = PairPotEiSlater1sp1spCorr(slater1s_widths, q2, q3, slater1p_widths, d2, d3, rcut, tr=tr)
    part_pair_slatercorr = ForcePartPair(system, nlist, scalings, pair_pot_slatercorr)

    # Outer atom loop
    for a0 in xrange(system.natom):
        # Inner atom loop; because of symmetry we only consider a1>a0.
        # The self-interactions for the case a1=a0 are handled later on.
        for a1 in xrange(a0+1,system.natom):
            #
            # Construct variable-variable tensor for this atom pair
            #
            tvv = np.zeros((nbv,nbv))
            for i in xrange(3):
                for j in xrange(3):
                    # Slater dipole - Slater dipole
                    pair_pot_point.dipoles[a0,i] = 1.0
                    pair_pot_slatercorr.slater1p_N[a0,i] = 1.0
                    pair_pot_point.dipoles[a1,j] = 1.0
                    pair_pot_slatercorr.slater1p_N[a1,j] = 1.0
                    tvv[i,j] = part_pair_point.compute() + part_pair_slatercorr.compute()
                    pair_pot_point.dipoles[:]  *=0.0
                    pair_pot_slatercorr.slater1p_N[:] *= 0.0
                    pair_pot_slatercorr.slater1s_N[:] *= 0.0
            # Fill big tensor
            Tvv[nbv*a0:nbv*(a0+1),nbv*a1:nbv*(a1+1)] = tvv.copy()
            Tvv[nbv*a1:nbv*(a1+1),nbv*a0:nbv*(a0+1)] = tvv.transpose().copy()
    return Tvv


def construct_ei_tensor_monodip(iterative):
    '''
    Construct the electrostatic interaction tensor with as basis functions:
        * a point monopole and a Slater monopole on the left
        * Slater dipoles with x, y and z orientation on the right
    '''
    system = iterative.ff.system
    nbf = 2
    nbv = 3
    Tfv = np.zeros((system.natom*nbf,system.natom*nbv)) # Tensor between fixed multipoles
    # Make some fake multipoles that help in calculating matrix elements
    # of the interaction tensor
    q1 = np.zeros((system.natom,))
    q2 = np.zeros((system.natom,))
    q3 = np.zeros((system.natom,))
    d1 = np.zeros((system.natom,3))
    d2 = np.zeros((system.natom,3))
    d3 = np.zeros((system.natom,3))
    # Construct potentials contributing to the tensor
    # Recycle settings from the force field
    alpha    = iterative.ff.part_pair_eidip.pair_pot.alpha
    rcut     = iterative.ff.part_pair_eidip.pair_pot.rcut
    tr       = iterative.ff.part_pair_eidip.pair_pot.get_truncation()
    nlist    = iterative.ff.part_pair_eidip.nlist
    scalings = iterative.ff.part_pair_eidip.scalings
    pair_pot_point  = PairPotEIDip(q1, d1, alpha, rcut, tr=tr )
    part_pair_point = ForcePartPair(system,  nlist, scalings, pair_pot_point)

    slater1s_widths = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_widths.copy()
    slater1p_widths = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1p_widths.copy()
    rcut            = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.rcut
    tr              = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.get_truncation()
    nlist           = iterative.ff.part_pair_eislater1sp1spcorr.nlist
    scalings        = iterative.ff.part_pair_eislater1sp1spcorr.scalings
    pair_pot_slatercorr  = PairPotEiSlater1sp1spCorr(slater1s_widths, q2, q3, slater1p_widths, d2, d3, rcut, tr=tr)
    part_pair_slatercorr = ForcePartPair(system, nlist, scalings, pair_pot_slatercorr)

    #
    # Construct fixed-variable tensor; not symmetric so apart from the rest
    #

    # Outer atom loop
    for a0 in xrange(system.natom):
        # Inner atom loop
        for a1 in xrange(system.natom):
            if a0==a1: continue
            #
            # Construct fixed-variable tensor for this atom pair
            #
            tfv = np.zeros((nbf,nbv))
            # Point charge - Slater dipole
            pair_pot_point.charges[a0] = 1.0
            pair_pot_slatercorr.slater1s_Z[a0] = 1.0
            for i in xrange(3):
                pair_pot_point.dipoles[a1,i] = 1.0
                pair_pot_slatercorr.slater1p_N[a1,i] = 1.0
                tfv[0,i] = part_pair_point.compute() + part_pair_slatercorr.compute()
                pair_pot_point.dipoles[:] *= 0.0
                pair_pot_slatercorr.slater1p_N[:] *= 0.0
            pair_pot_point.charges[:]  *= 0.0
            pair_pot_slatercorr.slater1s_Z[:] *= 0.0
            # Slater charge - Slater dipole
            pair_pot_point.charges[a0] = 1.0
            pair_pot_slatercorr.slater1s_N[a0] = 1.0
            for i in xrange(3):
                pair_pot_point.dipoles[a1,i] = 1.0
                pair_pot_slatercorr.slater1p_N[a1,i] = 1.0
                tfv[1,i] = part_pair_point.compute() + part_pair_slatercorr.compute()
                pair_pot_point.dipoles[:] *= 0.0
                pair_pot_slatercorr.slater1p_N[:] *= 0.0
            pair_pot_point.charges[:]  *= 0.0
            pair_pot_slatercorr.slater1s_N[:] *= 0.0                    
            # Fill big tensor
            Tfv[nbf*a0:nbf*(a0+1),nbv*a1:nbv*(a1+1)] = tfv.copy()
    return Tfv

class Polarization_FSPM_ISD(Hook):
    '''
    Fixed Slater and Point Monopole with Induced Slater Dipole. 
    '''
    def __init__(self, start=0, step=1, exclude=None):
        """
           **Arguments:**

           **Optional arguments:**

           start
                The first iteration at which this hook should be called.

           step
                The hook will be called every `step` iterations.

           exclude
                Specify list of pairs of atoms that should not induce each other.
                This can be used for instance to discard intramonomer polarization,
                by specifying a list of all atoms that are somehow connected.
        """
        self.exclude = exclude
        Hook.__init__(self, start, step)

    def __call__(self, iterative):
        # TODO: include periodic images in interaction tensor
        # TODO: this quickly becomes tedious because the computation of the
        # interaction tensor is not really compatible with the implementation
        # of the pair potentials.
        # The tensors should be computed by interfacing the low-level C energy
        # expression directly
        # TODO: what about contributions to the Hessian

        # Check that the necessary pair potentials are present in the force field.
        part_names = [part.name for part in iterative.ff.parts]
        assert 'pair_eidip' in part_names, "Forcefield has to contain pair_eidip when using FSPM_ISD polarizable model."
        assert 'pair_eislater1sp1spcorr' in part_names, "Forcefield has to contain pair_eislater1sp1spcorr when using FSPM_ISD polarizable model."
        assert 'polarization' in part_names, "Forcefield has to contain polarization term when using FSPM_ISD polarizable model."

        system = iterative.ff.system
        nbf = 2 # Number of basis functions for the fixed multipoles
        nbv = 3 # Number of basis functions for the induced multipoles

        #
        # Construct electrostatic interaction tensors
        #

        Tff = construct_ei_tensor_monomono(iterative)
        Tvv = construct_ei_tensor_dipdip(iterative)
        Tfv = construct_ei_tensor_monodip(iterative)

        #
        # Throw away unwanted interactions
        #

        if self.exclude is not None:
            for iatom0, iatom1 in self.exclude:
                Tff[nbf*iatom0:nbf*(iatom0+1),nbf*iatom1:nbf*(iatom1+1)] *= 0.0
                Tff[nbf*iatom1:nbf*(iatom1+1),nbf*iatom0:nbf*(iatom0+1)] *= 0.0
                Tvv[nbv*iatom0:nbv*(iatom0+1),nbv*iatom1:nbv*(iatom1+1)] *= 0.0
                Tvv[nbv*iatom1:nbv*(iatom1+1),nbv*iatom0:nbv*(iatom0+1)] *= 0.0
                Tfv[nbf*iatom0:nbf*(iatom0+1),nbv*iatom1:nbv*(iatom1+1)] *= 0.0
                Tfv[nbf*iatom1:nbf*(iatom1+1),nbv*iatom0:nbv*(iatom0+1)] *= 0.0

        #
        # Fill up diagonal elements of tensors
        #

        # The Tff diagonal elements do not matter; they are geometry independent
        # and only represent a constant offset of the total energy. WARNING:
        # this only holds if basis functions do not change during simulation
        # in our case this means the Slater widths need to be fixed.

        # The Tfv diagonal elements are all zero in our case, because the fixed
        # multipoles (and the their electrostatic potentials) are proportional
        # to Y_00 while the variable multipoles (and their potentials) are 
        # proportional to Y_1m. If both are on the same site, the integral
        # of rho*phi is zero.

        # The Tvv diagonal elements are the self-interaction of a Slater dipole
        # It is given by 7/32/24/width^3 TODO: check this expression
        for a0 in xrange(system.natom):
            for i in xrange(3):
                width = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1p_widths[a0,i]
                Tvv[nbv*a0+i,nbv*a0+i] = 7.0/32.0/24.0/width**3

        #
        # Calculate induced multipoles by solving linear system of equations
        #

        # Setup vector with fixed multipoles
        Qf = np.zeros((nbf*system.natom,))
        Qf[::2] = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_Z.copy()
        Qf[1::2] = iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1s_N.copy()

        # Simply invert the Tvv matrix
        #Qv = -np.dot(np.linalg.inv(Tvv), np.dot( Qf.transpose(), Tfv))
        # Solve iteratively
        Qv = np.zeros((nbv*system.natom,))
        converged = False
        maxiterations = 100
        threshold = 1e-6
        niter = 0
        while not converged and niter<maxiterations:
            Qv_new = -np.dot(Qf.transpose(),Tfv)
            for i in xrange(nbv*system.natom):
                for j in xrange(nbv*system.natom):
                    if i==j: continue
                    Qv_new[i] -= Tvv[i,j]*Qv[j]
                Qv_new[i] /= Tvv[i,i]
            difference = np.amax(Qv-Qv_new)
            if difference < threshold: converged = True
            Qv = Qv_new.copy()
            niter += 1

        #
        # Calculate different parts of the total energy
        #

        # Calculate energy it takes to polarize the dipoles
        Epol = 0.0
        for a0 in xrange(system.natom):
            for i in xrange(3):
                Epol += 0.5*Tvv[nbv*a0+i,nbv*a0+i]*Qv[nbv*a0+i]**2
        Eff = 0.5*np.dot(np.dot(Qf.transpose(),Tff),Qf)
        Evv = 0.5*np.dot(np.dot(Qv.transpose(),Tvv),Qv) - Epol
        Efv = 1.0*np.dot(np.dot(Qf.transpose(),Tfv),Qv)
        # This is the complete electrostatic interaction, including self-interactions
        # (except for the fixed multipoles) 
        Etot = Eff + Evv + Efv + Epol
        # This is the difference in electrostatic interaction energy (excluding
        # self-interactions) by inducing multipoles.
        Eind = Evv + Efv
        # Check relationship between different parts of the total energy
        assert np.abs(Eind+2.0*Epol+Evv) < 1e-4
        
        #
        # Update dipoles
        #
        iterative.ff.part_pair_eislater1sp1spcorr.pair_pot.slater1p_N[:] = Qv.reshape((system.natom,3))
        iterative.ff.part_pair_eidip.pair_pot.dipoles[:] = Qv.reshape((system.natom,3))
        
        # Add polarization energy to the force field
        iterative.ff.part_polarization.set_pol_energy(Epol)
