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

from molmod import check_delta

from yaff import *


__all__ = [
    'check_gpos_part', 'check_vtens_part', 'check_hess_part',
    'check_gpos_ff', 'check_vtens_ff',
]


def check_gpos_part(system, part, nlists=None):
    def fn(x, do_gradient=False):
        system.pos[:] = x.reshape(system.natom, 3)
        if nlists is not None:
            nlists.update()
        if do_gradient:
            gpos = np.zeros(system.pos.shape, float)
            e = part.compute(gpos)
            assert np.isfinite(e)
            assert np.isfinite(gpos).all()
            return e, gpos.ravel()
        else:
            e = part.compute()
            assert np.isfinite(e)
            return e

    x = system.pos.ravel()
    dxs = np.random.normal(0, 1e-4, (100, len(x)))
    check_delta(fn, x, dxs)


def check_vtens_part(system, part, nlists=None, symm_vtens=True):
    '''
        * symm_vtens: Perform check that the virial tensor is a symmetric matrix.
                      For instance for dipole interactions, this can be not the case.
    '''
    # define some rvecs and gvecs
    if system.cell.nvec == 3:
        gvecs = system.cell.gvecs
        rvecs = system.cell.rvecs
    else:
        gvecs = np.identity(3, float)
        rvecs = np.identity(3, float)

    # Get the reduced coordinates
    reduced = np.dot(system.pos, gvecs.transpose())
    if symm_vtens: assert abs(np.dot(reduced, rvecs) - system.pos).max() < 1e-10

    def fn(x, do_gradient=False):
        rvecs = x.reshape(3, 3)
        if system.cell.nvec == 3:
            system.cell.update_rvecs(rvecs)
        system.pos[:] = np.dot(reduced, rvecs)
        if nlists is not None:
            nlists.update()
        if do_gradient:
            vtens = np.zeros((3, 3), float)
            e = part.compute(vtens=vtens)
            gvecs = np.linalg.inv(rvecs).transpose()
            grvecs = np.dot(gvecs, vtens)
            assert np.isfinite(e)
            assert np.isfinite(vtens).all()
            assert np.isfinite(grvecs).all()
            if symm_vtens: assert abs(vtens - vtens.transpose()).max() < 1e-10
            return e, grvecs.ravel()
        else:
            e = part.compute()
            assert np.isfinite(e)
            return e

    x = rvecs.ravel()
    dxs = np.random.normal(0, 1e-4, (100, len(x)))
    check_delta(fn, x, dxs)


def check_hess_part(system, part, nlists=None, numgrad=False):
    def fn(x, do_gradient=False):
        system.pos[:] = x.reshape(system.natom, 3)
        if nlists is not None:
            nlists.update()
        if do_gradient:
            gpos = np.zeros(system.pos.shape, float)
            hess = np.zeros((np.prod(gpos.shape),np.prod(gpos.shape)),float)
            e = part.compute(gpos,hess=hess)
            assert np.isfinite(e)
            assert np.isfinite(gpos).all()
            assert np.isfinite(hess).all()
            return gpos.ravel()[iatom], hess[iatom,:]
        else:
            gpos = np.zeros(system.pos.shape, float)
            e = part.compute(gpos)
            assert np.isfinite(e)
            assert np.isfinite(gpos).all()
            return gpos.ravel()[iatom]

    if numgrad:
        x_copy = system.pos.copy()
        x = system.pos.ravel()
        dxs = np.random.normal(0, 1e-4, (100, len(x)))
        # Column i of the hessian is the derivative of element i of the gradient.
        # This is checked by computing that derivative numerically
        # This test can take a while...
        for iatom in xrange(3*system.natom):
            check_delta(fn, x, dxs)
        system.pos[:] = x_copy
    # Next we compare the analytical hessian with the numerical estimate from
    # Yaff. This is actually more or less the same as above...
    # Values of the thresholds are debatable...
    if nlists is not None:
        nlists.update()
    hessian_ana = np.zeros((np.prod(system.pos.shape),np.prod(system.pos.shape)),float)
    part.compute(hess=hessian_ana)
    hessian_num = estimate_cart_hessian(ForceField(system, [part], nlists))
    print hessian_ana[:6,:6]
    print hessian_num[:6,:6]
    av_err = np.mean((hessian_num-hessian_ana)**2)
    # To compute a meaningful relative error, mask out very small reference values
    mask = np.abs(hessian_num) > 1e-8
    max_relerr = np.amax( (hessian_ana[mask]/hessian_num[mask]-1.0)**2)
    index_relerr = np.argmax( (hessian_ana[mask]/hessian_num[mask]-1.0)**2)
    assert av_err<1e-10, "Average hessian error too large %5.1e" % av_err
    assert max_relerr<1e-7, "Largest relative error of hessian too large %5.1e (%+5.1e - %+5.1e)" % (max_relerr,hessian_num[mask][index_relerr],hessian_ana[mask][index_relerr])


def check_gpos_ff(ff):
    def fn(x, do_gradient=False):
        ff.update_pos(x.reshape(ff.system.natom, 3))
        if do_gradient:
            gpos = np.zeros(ff.system.pos.shape, float)
            e = ff.compute(gpos)
            assert np.isfinite(e)
            assert np.isfinite(gpos).all()
            return e, gpos.ravel()
        else:
            e = ff.compute()
            assert np.isfinite(e)
            return e

    x = ff.system.pos.ravel()
    dxs = np.random.normal(0, 1e-4, (100, len(x)))
    check_delta(fn, x, dxs)


def check_vtens_ff(ff):
    # define some rvecs and gvecs
    if ff.system.cell.nvec == 3:
        gvecs = ff.system.cell.gvecs
        rvecs = ff.system.cell.rvecs
    else:
        gvecs = np.identity(3, float)
        rvecs = np.identity(3, float)

    # Get the reduced coordinates
    reduced = np.dot(ff.system.pos, gvecs.transpose())
    assert abs(np.dot(reduced, rvecs) - ff.system.pos).max() < 1e-10

    def fn(x, do_gradient=False):
        rvecs = x.reshape(3, 3)
        if ff.system.cell.nvec == 3:
            ff.update_rvecs(rvecs)
        ff.update_pos(np.dot(reduced, rvecs))
        if do_gradient:
            vtens = np.zeros((3, 3), float)
            e = ff.compute(vtens=vtens)
            gvecs = np.linalg.inv(rvecs).transpose()
            grvecs = np.dot(gvecs, vtens)
            assert np.isfinite(e)
            assert np.isfinite(vtens).all()
            assert np.isfinite(grvecs).all()
            assert abs(vtens - vtens.transpose()).max() < 1e-10
            return e, grvecs.ravel()
        else:
            e = ff.compute()
            assert np.isfinite(e)
            return e

    x = rvecs.ravel()
    dxs = np.random.normal(0, 1e-4, (100, len(x)))
    check_delta(fn, x, dxs)
