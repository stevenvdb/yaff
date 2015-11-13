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
//#define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif

#include <cmath>
#include <cstring>
#include <stdexcept>

#include <iostream>
#include <iomanip>
#include "pair_pot.h"
#include "nlist.h"
#include "splinepot.h"
#include <omp.h>
#include <map>


SplinePot::SplinePot(std::size_t n, long *ffatype_ids, long nffa, double *r, double *phi, double *dphi, long npoints, scaling_row_type *stab,
                        long nstab, long natom):
    r(r), phi(phi), dphi(dphi), npoints(npoints), nffa(nffa), ffatype_ids(ffatype_ids)
{
    first_x = r[0];
    last_x = r[npoints-1];
    alpha = (last_x-first_x)/(npoints-1);
    // Construct map containing scaling of 1-2,1-3 and 1-4 neighbors
    long i;
    for (i=0; i<nstab; i++) {
      //printf("%5d %5d %8d %6.3f\n", stab[i].a, stab[i].b,stab[i].b+ nneigh*stab[i].a, stab[i].scale);
      mymap[stab[i].b + natom*stab[i].a] = stab[i].scale;
      mymap[stab[i].b*natom + stab[i].a] = stab[i].scale;
  }
}


double SplinePot::eval(double x, long index ) {
    // 1) find the index of the interval in which t lies.
    double t = (x-first_x)/alpha;
    int j = (int)floor(t);
    if (j==npoints - 1) j = npoints - 2;
    // 2) do the interpolation
    double u = t - j;
    j += npoints*index;
    double z = phi[j+1] - phi[j];
    return phi[j] + u*(dphi[j] + u*(3*z - 2*dphi[j] - dphi[j+1] + u*(-2*z + dphi[j] + dphi[j+1])));
}


double SplinePot::eval_deriv(double x, long index ) {
    // 1) find the index of the interval in which t lies.
    double t = (x-first_x)/alpha;
    int j = (int)floor(t);
    if (j==npoints - 1) j = npoints - 2;
    // 2) do the interpolation
    double u = t - j;
    j += npoints*index;
    double z = phi[j+1] - phi[j];
    return (dphi[j] + u*(6*z - 4*dphi[j] - 2*dphi[j+1] + u*(-6*z + 3*dphi[j] + 3*dphi[j+1])))/alpha;
}


double SplinePot::spline_pot_compute(neigh_row_type *neighs,
                        long nneigh, double *gpos, double* vtens, long ndof, long natom)
{
  long i, center_index, other_index, index, iffa, jffa;
  double s, energy, vg, h;
  energy = 0.0;
  // Create iterator that will be used to determine scaling of interactions
  std::map<long, double>::iterator it;
  #pragma omp parallel private(i)
  {
    // Construct a private copy of gpos for every thread
    double gpos_private[ndof];
    for (i=0; i<ndof; i++) gpos_private[i] = 0.0;
    double vtens_private[9];
    for (i=0; i<9; i++) vtens_private[i] = 0.0;
    // Compute the interactions.
    #pragma omp for reduction(+:energy) ordered private(center_index,other_index, index, iffa, jffa, s, h, vg) firstprivate(it) schedule(static)
    for (i=0; i<nneigh; i++) {
      // Find the scale
      center_index = neighs[i].a;
      other_index = neighs[i].b;
      if ((neighs[i].r0 == 0) && (neighs[i].r1 == 0) && (neighs[i].r2 == 0)) {
        it = mymap.find(natom*center_index+other_index);
        if (it == mymap.end()) s = 1.0;
        else s = it->second;
      } else {
        s = 1.0;
      }
      // If the scale is non-zero, compute the contribution.
      if (s > 0.0) {
        // Find the index of the spline representing the interaction between
        // these atom types
        iffa = ffatype_ids[center_index];
        jffa = ffatype_ids[other_index];
        if ( iffa > jffa ) {
          jffa = ffatype_ids[center_index];
          iffa = ffatype_ids[other_index];
        }
        index = nffa*iffa - (iffa*(iffa-1))/2 + jffa - iffa;
        // Evaluate the spline at the current interatomic distance
        energy += s*eval(neighs[i].d, index);
        if ((gpos!=NULL) || (vtens!=NULL)) vg = s*eval_deriv(neighs[i].d, index)/neighs[i].d;
        if (gpos!=NULL) {
          h = neighs[i].dx*vg;
          gpos_private[3*other_index  ] += h;
          gpos_private[3*center_index   ] -= h;
          h = neighs[i].dy*vg;
          gpos_private[3*other_index+1] += h;
          gpos_private[3*center_index +1] -= h;
          h = neighs[i].dz*vg;
          gpos_private[3*other_index+2] += h;
          gpos_private[3*center_index +2] -= h;
        }
        if (vtens!=NULL) {
            vtens_private[0] += neighs[i].dx*neighs[i].dx*vg;
            vtens_private[4] += neighs[i].dy*neighs[i].dy*vg;
            vtens_private[8] += neighs[i].dz*neighs[i].dz*vg;
            h = neighs[i].dx*neighs[i].dy*vg;
            vtens_private[1] += h;
            vtens_private[3] += h;
            h = neighs[i].dx*neighs[i].dz*vg;
            vtens_private[2] += h;
            vtens_private[6] += h;
            h = neighs[i].dy*neighs[i].dz*vg;
            vtens_private[5] += h;
            vtens_private[7] += h;
        }
      }
    }
    #pragma omp critical
    {
      if (gpos!=NULL) {
        for(i=0; i<ndof; i++)
        {
          gpos[i] += gpos_private[i];
        }
      }
      if (vtens!=NULL) {
        for(i=0; i<9; i++)
        { vtens[i] += vtens_private[i];
        }
      }
    }
  }
  return energy;
}
