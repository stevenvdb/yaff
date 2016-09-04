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
    'check_gpos_ff', 'check_vtens_ff', 'check_hess_ff', 'check_jacobian_ic',
    'check_hessian_oneic',
]


def check_gpos_part(system, part, nlists=None):
    def fn(x, do_gradient=False):
        system.pos[:] = x.reshape(system.natom, 3)
        if nlists is not None:
            nlists.update()
        if do_gradient:
            gpos = np.zeros(system.pos.shape, float)
            hess = np.zeros((3*system.natom,3*system.natom), float)
            e = part.compute(gpos, hess=hess)
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


def check_hess_numgrad_part(system, part, nlists=None):
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

    x = system.pos.ravel()
    dxs = np.random.normal(0, 1e-4, (100, len(x)))
    # Column i of the hessian is the derivative of element i of the gradient.
    # This is checked by computing that derivative numerically
    # This test can take a while...
    for iatom in xrange(3*system.natom):
        check_delta(fn, x, dxs)


def check_hess_part(system, part, nlists=None):
    # We compare the analytical hessian with the numerical estimate from
    # Yaff. This is actually more or less the same as above...
    # Values of the thresholds are debatable...
    if nlists is not None:
        nlists.update()
    hessian_ana = np.zeros((np.prod(system.pos.shape),np.prod(system.pos.shape)),float)
    gpos = np.zeros(system.pos.shape, float)
    part.compute(hess=hessian_ana)
    hessian_num = estimate_cart_hessian(ForceField(system, [part], nlists),eps=1e-6)
    check_hess(hessian_ana, hessian_num)


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


def check_hess(hessian_ana, hessian_num):
    av_err = np.mean((hessian_num-hessian_ana)**2)
    # To compute a meaningful relative error, mask out very small reference values
    mask = np.abs(hessian_num) > 1e-7
    #    for i in xrange(hessian_ana.shape[0]):
    #        for j in xrange(hessian_ana.shape[0]):
    #            print "%3d %3d %+12.2e %+12.2e" % (i,j,hessian_ana[i,j],hessian_num[i,j])
    if np.sum(mask)>0:
        max_relerr = np.amax( (hessian_ana[mask]/hessian_num[mask]-1.0)**2)
        index_relerr = np.argmax( (hessian_ana[mask]/hessian_num[mask]-1.0)**2)
        assert av_err<1e-10, "Average hessian error too large %5.1e" % av_err
        assert max_relerr<1e-6, "Largest relative error of hessian too large %5.1e (%+5.1e - %+5.1e)" % (max_relerr,hessian_num[mask][index_relerr],hessian_ana[mask][index_relerr])


def check_hess_ff(ff):
    hessian_ana = np.zeros((np.prod(ff.system.pos.shape),np.prod(ff.system.pos.shape)),float)
    ff.compute(hess=hessian_ana)
    hessian_num = estimate_cart_hessian(ff)
    check_hess(hessian_ana, hessian_num)


def check_jacobian_ic(system, dlist, iclist, iic=0):
    '''
    Check the implementation of the derivative of one internal coordinate
    towards cartesian coordinates by comparison with the result obtained using
    finite differences.
    '''
    def fn(x, do_gradient=False):
        system.pos[:] = x.reshape(system.natom, 3)
        dlist.forward()
        iclist.forward()
        if do_gradient:
            q = iclist.ictab[iic]['value']
            jacobian = np.zeros(np.prod(system.pos.shape))
            ic_jac = iclist.jacobian()[iic]
            kind = iclist.ictab[iic]['kind']
            # Loop over all relative vectors making up this internal coordinate
            index = iclist.ictab[iic]['i0']
            i,j = dlist.deltas[index]['i'],dlist.deltas[index]['j']
            d = iclist.ictab[iic]['i0']
            jacobian[3*i:3*(i+1)] -= iclist.ictab[iic]['sign0']**0*ic_jac[3*d:3*(d+1)]
            jacobian[3*j:3*(j+1)] += iclist.ictab[iic]['sign0']**0*ic_jac[3*d:3*(d+1)]
            if not kind in [0,5]:
                index = iclist.ictab[iic]['i1']
                i,j = dlist.deltas[index]['i'],dlist.deltas[index]['j']
                d = iclist.ictab[iic]['i1']
                jacobian[3*i:3*(i+1)] -= iclist.ictab[iic]['sign1']**0*ic_jac[3*d:3*(d+1)]
                jacobian[3*j:3*(j+1)] += iclist.ictab[iic]['sign1']**0*ic_jac[3*d:3*(d+1)]
            if not kind in [0,1,2,5]:
                index = iclist.ictab[iic]['i2']
                i,j = dlist.deltas[index]['i'],dlist.deltas[index]['j']
                d = iclist.ictab[iic]['i2']
                jacobian[3*i:3*(i+1)] -= iclist.ictab[iic]['sign2']**0*ic_jac[3*d:3*(d+1)]
                jacobian[3*j:3*(j+1)] += iclist.ictab[iic]['sign2']**0*ic_jac[3*d:3*(d+1)]
            return q, jacobian
        else:
            return iclist.ictab[iic]['value']
    x = system.pos.ravel()
    dxs = np.random.normal(0, 1e-4, (100, len(x)))
    check_delta(fn, x, dxs)


def check_hessian_oneic(dlist, iclist):
    '''Check the second derivate of an internal coordinate towards the delta
    vectors involved for this internal coordinate by comparing with finite
    difference approach
    '''
    class OneIC(object):
        def __init__(self, dlist, iclist):
            assert iclist.nic == 1
            self.dlist = dlist
            self.iclist = iclist
            self.dlist.forward()
            self.x0 = np.zeros((dlist.ndelta*3,))
            for i in xrange(dlist.ndelta):
                self.x0[3*i+0] = self.dlist.deltas[i]['dx'].copy()
                self.x0[3*i+1] = self.dlist.deltas[i]['dy'].copy()
                self.x0[3*i+2] = self.dlist.deltas[i]['dz'].copy()

        def fun(self, x, do_gradient=False):
            for i in xrange(dlist.ndelta):
                self.dlist.deltas[i]['dx'] = x[3*i+0]
                self.dlist.deltas[i]['dy'] = x[3*i+1]
                self.dlist.deltas[i]['dz'] = x[3*i+2]
            self.iclist.forward()
            q = self.iclist.ictab[0]['value']
            gradient = self.iclist.jacobian()
            return q, gradient

        def reset(self):
            for i in xrange(dlist.ndelta):
                self.dlist.deltas[i]['dx'] = self.x0[3*i+0]
                self.dlist.deltas[i]['dy'] = self.x0[3*i+1]
                self.dlist.deltas[i]['dz'] = self.x0[3*i+2]

    def estimate_hessian(oneic, eps=1e-4):
        # Loop over all displacements
        x1 = oneic.x0.copy()
        rows = np.zeros((len(x1), len(x1)), float)
        for i in xrange(len(x1)):
            x1[i] = oneic.x0[i] + eps
            epot, gradient_p = oneic.fun(x1, do_gradient=True)
            x1[i] = oneic.x0[i] - eps
            epot, gradient_m = oneic.fun(x1, do_gradient=True)
            rows[i] = (gradient_p-gradient_m)/(2*eps)
            x1[i] = oneic.x0[i]
        oneic.reset()
        # Enforce symmetry and return
        return 0.5*(rows + rows.T)

    # Estimate hessian with finite difference
    oneic = OneIC(dlist, iclist)
    hessian_num = estimate_hessian(oneic)
    oneic.reset()
    # Analytical hessian
    dlist.forward()
    iclist.forward()
    iclist.ictab[0]['grad'] = 1.0
    hessian_ana = iclist.hessian()
    check_hess(hessian_ana, hessian_num)
