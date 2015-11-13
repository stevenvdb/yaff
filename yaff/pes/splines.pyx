#!/usr/bin/env python
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
'''C++ extension'''

import numpy as np

cimport numpy as np

cimport nlist
cimport pair_pot
cimport splines

__all__ = ['PairPotSpline']

cdef class PairPotSpline:
    '''A tabulated pair potential with interpolation using cubic splines
       supporting OpenMP parallelism.

       **Arguments:**

       ffatype_ids
            NumPy array containing the atom-type identifiers

       nffa
            The number of different atom types

       r
            NumPy array of equidistant radial grid points

       phi
            NumPy array of the potential evaluated at the radial grid points

       dphi
            NumPy array of the derivative of the potential evaluated at the
            radial grid points
    '''
    cdef splines.SplinePot* _this
    cdef np.ndarray _r
    cdef np.ndarray _phi
    cdef np.ndarray _dphi
    cdef long _natom

    def __cinit__(self, long nsplines,
                np.ndarray[long, ndim=1] ffatype_ids, long nffa,
                np.ndarray[double, ndim=1] r not None,
                np.ndarray[double, ndim=2] phi not None,
                np.ndarray[double, ndim=2] dphi not None,
                np.ndarray[pair_pot.scaling_row_type, ndim=1] stab):
        assert ffatype_ids.flags['C_CONTIGUOUS']
        assert r.flags['C_CONTIGUOUS']
        assert phi.flags['C_CONTIGUOUS']
        assert dphi.flags['C_CONTIGUOUS']

        self._r = r
        self._phi = phi
        self._dphi = dphi
        self._natom = len(ffatype_ids)

        # Check that r is equidistantly spaced
        diff = r[1] - r[0]
        assert np.all( np.abs( r[1:]-r[:-1]-diff ) < 1e-10), "Data points of r are not equidistant"
        # Number of grid points
        n = r.shape[0]
        # Check that phi and dphi contain proper amount of points
        assert phi.shape[1]==n
        assert dphi.shape[1]==n
        # Rescale gradient
        self._dphi *= (r[-1]-r[0])/(n-1)
        # Construct C++ object
        self._this = new splines.SplinePot(nsplines,  <long*>ffatype_ids.data, nffa, <double*>r.data,
                <double*>phi.data, <double*>dphi.data, n, <pair_pot.scaling_row_type*>stab.data, len(stab), len(ffatype_ids))

    def __dealloc__(self):
        del self._this

    def compute(self, np.ndarray[nlist.neigh_row_type, ndim=1] neighs,
                np.ndarray[double, ndim=2] gpos,
                np.ndarray[double, ndim=2] vtens, long nneigh):
        cdef double *my_gpos
        cdef double *my_vtens

        assert neighs.flags['C_CONTIGUOUS']

        if gpos is None:
            my_gpos = NULL
            ndof = 0
        else:
            assert gpos.flags['C_CONTIGUOUS']
            assert gpos.shape[1] == 3
            my_gpos = <double*>gpos.data
            ndof = gpos.shape[0]*gpos.shape[1]

        if vtens is None:
            my_vtens = NULL
        else:
            assert vtens.flags['C_CONTIGUOUS']
            assert vtens.shape[0] == 3
            assert vtens.shape[1] == 3
            my_vtens = <double*>vtens.data

        return self._this.spline_pot_compute(<nlist.neigh_row_type*>neighs.data, nneigh,
                my_gpos, my_vtens, ndof, self._natom
            )
