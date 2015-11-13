// YAFF is yet another force-field code
// Copyright (C) 2011 - 2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
// Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>, Center for Molecular Modeling
// (CMM), Ghent University, Ghent, Belgium; all rights reserved unless otherwise
// stated.
//
// This file is part of YAFF.
//
// YAFF is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// YAFF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
//--
#include "pair_pot.h"
#include "nlist.h"
#include <vector>
#include <map>

#ifndef SPLINEPOT_H
#define SPLINEPOT_H

class SplinePot {
    private:
        double * r;
        double * phi;
        double * dphi;
        long npoints;
        double first_x;
        double last_x;
        double alpha;
        std::map<long, double> mymap;
    public:
        long nffa;
        long *ffatype_ids;
        explicit SplinePot(std::size_t n, long *ffatype_ids, long nffa, double *r, double *phi, double *dphi, long npoints, scaling_row_type *stab,
                        long nstab, long natom);
        double eval(double x, long index);
        double eval_deriv(double x, long index);
        double spline_pot_compute (neigh_row_type *neighs,
                        long nneigh, double *gpos, double* vtens, long ndof, long natom);
};

#endif
