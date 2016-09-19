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

#include <math.h>
#include "iclist.h"
#include "constants.h"
#include <stdio.h>

typedef double (*ic_forward_type)(iclist_row_type*, dlist_row_type*);

double forward_bond(iclist_row_type* ic, dlist_row_type* deltas) {
  double *delta;
  delta = (double*)(deltas + (*ic).i0);
  return sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
}

double forward_bend_cos(iclist_row_type* ic, dlist_row_type* deltas) {
  double *delta0, *delta1;
  double d0, d1, dot;
  delta0 = (double*)(deltas + (*ic).i0);
  delta1 = (double*)(deltas + (*ic).i1);
  d0 = sqrt(delta0[0]*delta0[0] + delta0[1]*delta0[1] + delta0[2]*delta0[2]);
  d1 = sqrt(delta1[0]*delta1[0] + delta1[1]*delta1[1] + delta1[2]*delta1[2]);
  if ((d0 == 0) || (d1 == 0)) return 0.0;
  dot = delta0[0]*delta1[0] + delta0[1]*delta1[1] + delta0[2]*delta1[2];
  return (*ic).sign0*(*ic).sign1*dot/d0/d1;
}

double forward_bend_angle(iclist_row_type* ic, dlist_row_type* deltas) {
  double c;
  c = forward_bend_cos(ic, deltas);
  if (c>1.0) c = 1.0;
  if (c<-1.0) c = -1.0;
  return acos(c);
}

double forward_dihed_cos(iclist_row_type* ic, dlist_row_type* deltas) {
  long i;
  double *delta0, *delta1, *delta2;
  double a[3], b[3];
  double tmp0, tmp1, tmp2;
  delta0 = (double*)(deltas + (*ic).i0);
  delta1 = (double*)(deltas + (*ic).i1);
  delta2 = (double*)(deltas + (*ic).i2);
  tmp1 = sqrt(delta1[0]*delta1[0] + delta1[1]*delta1[1] + delta1[2]*delta1[2]);
  tmp0 = (delta0[0]*delta1[0] + delta0[1]*delta1[1] + delta0[2]*delta1[2])/tmp1;
  tmp2 = (delta1[0]*delta2[0] + delta1[1]*delta2[1] + delta1[2]*delta2[2])/tmp1;
  for (i=0; i<3; i++) {
    a[i] = delta0[i] - tmp0*delta1[i]/tmp1;
    b[i] = delta2[i] - tmp2*delta1[i]/tmp1;
  }
  tmp0 = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  tmp2 = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
  tmp1 = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return (*ic).sign0*(*ic).sign2*tmp1/tmp0/tmp2;
}

double forward_dihed_angle(iclist_row_type* ic, dlist_row_type* deltas) {
  double c;
  c = forward_dihed_cos(ic, deltas);
  // Guard against round-off errors before taking the dot product.
  if (c > 1) {
    c = 1;
  } else if (c < -1) {
    c = -1;
  }
  return acos(c);
}

double forward_oop_cos_low(double *delta0, double *delta1, double *delta2) {
  double n[3];
  double n_sq, tmp0, tmp1;
  // The normal to the plane spanned by the first and second vector
  n[0] = delta0[1]*delta1[2] - delta0[2]*delta1[1];
  n[1] = delta0[2]*delta1[0] - delta0[0]*delta1[2];
  n[2] = delta0[0]*delta1[1] - delta0[1]*delta1[0];
  // The norm squared of this normal
  n_sq = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  // Norm squared of the third vector
  tmp0 = delta2[0]*delta2[0] + delta2[1]*delta2[1] + delta2[2]*delta2[2];
  // Dot product of the normal with the third vector
  tmp1 = n[0]*delta2[0] + n[1]*delta2[1] + n[2]*delta2[2];
  // The cosine of the oop angle (assumed to be positive)
  return sqrt(1.0 - tmp1*tmp1/tmp0/n_sq);
}

double forward_oop_cos(iclist_row_type* ic, dlist_row_type* deltas) {
  double *delta0, *delta1, *delta2;
  delta0 = (double*)(deltas + (*ic).i0);
  delta1 = (double*)(deltas + (*ic).i1);
  delta2 = (double*)(deltas + (*ic).i2);
  return forward_oop_cos_low(delta0, delta1, delta2);
}

double forward_oop_meancos(iclist_row_type* ic, dlist_row_type* deltas) {
  double *delta0, *delta1, *delta2;
  double tmp;
  delta0 = (double*)(deltas + (*ic).i0);
  delta1 = (double*)(deltas + (*ic).i1);
  delta2 = (double*)(deltas + (*ic).i2);
  tmp = forward_oop_cos_low(delta0, delta1, delta2);
  tmp += forward_oop_cos_low(delta2, delta0, delta1);
  tmp += forward_oop_cos_low(delta1, delta2, delta0);
  return tmp/3;
}

double forward_oop_angle_low(double *delta0, double *delta1, double *delta2) {
  double c;
  c = forward_oop_cos_low(delta0, delta1, delta2);
  // Guard against round-off errors before taking the dot product.
  if (c > 1) {
    c = 1;
  } else if (c < -1) {
    c = -1;
  }
  return acos(c);
}

double forward_oop_angle(iclist_row_type* ic, dlist_row_type* deltas) {
  double *delta0, *delta1, *delta2;
  delta0 = (double*)(deltas + (*ic).i0);
  delta1 = (double*)(deltas + (*ic).i1);
  delta2 = (double*)(deltas + (*ic).i2);
  return forward_oop_angle_low(delta0, delta1, delta2);
}

double forward_oop_meanangle(iclist_row_type* ic, dlist_row_type* deltas) {
  double *delta0, *delta1, *delta2;
  double tmp;
  delta0 = (double*)(deltas + (*ic).i0);
  delta1 = (double*)(deltas + (*ic).i1);
  delta2 = (double*)(deltas + (*ic).i2);
  tmp = forward_oop_angle_low(delta0, delta1, delta2);
  tmp += forward_oop_angle_low(delta2, delta0, delta1);
  tmp += forward_oop_angle_low(delta1, delta2, delta0);
  return tmp/3;
}

double forward_oop_distance(iclist_row_type* ic, dlist_row_type* deltas) {
  double *delta0, *delta1, *delta2;
  double n[3];
  double n_norm, n_dot_d2;
  delta0 = (double*)(deltas + (*ic).i0);
  delta1 = (double*)(deltas + (*ic).i1);
  delta2 = (double*)(deltas + (*ic).i2);
  // The normal to the plane spanned by the first and second vector
  n[0] = delta0[1]*delta1[2] - delta0[2]*delta1[1];
  n[1] = delta0[2]*delta1[0] - delta0[0]*delta1[2];
  n[2] = delta0[0]*delta1[1] - delta0[1]*delta1[0];
  // The norm of this normal
  n_norm = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );
  // If n_norm==0, the first and second vector don't span a plane!
  if (n_norm == 0) return 0.0;
  n_dot_d2 = n[0]*delta2[0] + n[1]*delta2[1] + n[2]*delta2[2];
  // Distance from point to plane spanned by first and second vector
  return n_dot_d2/n_norm*(*ic).sign0*(*ic).sign1*(*ic).sign2;
}

ic_forward_type ic_forward_fns[11] = {
  forward_bond, forward_bend_cos, forward_bend_angle, forward_dihed_cos, forward_dihed_angle, forward_bond,
  forward_oop_cos, forward_oop_meancos, forward_oop_angle, forward_oop_meanangle, forward_oop_distance
};

void iclist_forward(dlist_row_type* deltas, iclist_row_type* ictab, long nic) {
  long i;
  for (i=0; i<nic; i++) {
    ictab[i].value = ic_forward_fns[ictab[i].kind](ictab + i, deltas);
    ictab[i].grad = 0.0;
  }
}


typedef void (*ic_back_type)(iclist_row_type*, dlist_row_type*, double, double, double*);

void back_bond(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double *avtens) {
  dlist_row_type *delta;
  double x, fx0, fy0, fz0;
  delta = deltas + (*ic).i0;
  x = grad/value;
  fx0 = x*(*delta).dx;
  fy0 = x*(*delta).dy;
  fz0 = x*(*delta).dz;
  (*delta).gx += fx0;
  (*delta).gy += fy0;
  (*delta).gz += fz0;
  if (avtens != NULL) {
      // Forces from i on j, factor 0.5 for two-body
      fx0 *= 0.5*(*ic).sign0;
      fy0 *= 0.5*(*ic).sign0;
      fz0 *= 0.5*(*ic).sign0;
      // Figure out (i,j,k) indexes of particles present in ic
      long indexes[2], i;
      indexes[0] = (*delta).i;
      indexes[1] = (*delta).j;
      printf("INDEXES %5ld %5ld\n",indexes[0],indexes[1]);
      // Fill atomic virial tensors for every particle in ic
      for (i=0;i<2;i++) {
        avtens[indexes[i]*9  ] += fx0*(*delta).dx;
        avtens[indexes[i]*9+1] += fx0*(*delta).dy;
        avtens[indexes[i]*9+2] += fx0*(*delta).dz;
        avtens[indexes[i]*9+3] += fy0*(*delta).dx;
        avtens[indexes[i]*9+4] += fy0*(*delta).dy;
        avtens[indexes[i]*9+5] += fy0*(*delta).dz;
        avtens[indexes[i]*9+6] += fz0*(*delta).dx;
        avtens[indexes[i]*9+7] += fz0*(*delta).dy;
        avtens[indexes[i]*9+8] += fz0*(*delta).dz;
      }
  }
}

void back_bend_cos(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double *avtens) {
  dlist_row_type *delta0, *delta1;
  double e0[3], e1[3];
  double d0, d1, fac;
  double fx0, fx1, fy0, fy1, fz0, fz1;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  d0 = sqrt((*delta0).dx*(*delta0).dx + (*delta0).dy*(*delta0).dy + (*delta0).dz*(*delta0).dz);
  d1 = sqrt((*delta1).dx*(*delta1).dx + (*delta1).dy*(*delta1).dy + (*delta1).dz*(*delta1).dz);
  e0[0] = (*delta0).dx/d0;
  e0[1] = (*delta0).dy/d0;
  e0[2] = (*delta0).dz/d0;
  e1[0] = (*delta1).dx/d1;
  e1[1] = (*delta1).dy/d1;
  e1[2] = (*delta1).dz/d1;
  fac = (*ic).sign0*(*ic).sign1;
  grad *= fac;
  value *= fac;
  fac = grad/d0;
  fx0 = fac*(e1[0] - value*e0[0]);
  fy0 = fac*(e1[1] - value*e0[1]);
  fz0 = fac*(e1[2] - value*e0[2]);
  (*delta0).gx += fx0;
  (*delta0).gy += fy0;
  (*delta0).gz += fz0;
  fac = grad/d1;
  fx1 = fac*(e0[0] - value*e1[0]);
  fy1 = fac*(e0[1] - value*e1[1]);
  fz1 = fac*(e0[2] - value*e1[2]);
  (*delta1).gx += fx1;
  (*delta1).gy += fy1;
  (*delta1).gz += fz1;
  if (avtens != NULL) {
      // Forces from i on j, factor M_THIRD for three-body
      fx0 *= M_THIRD;//*(*ic).sign0;
      fy0 *= M_THIRD;//*(*ic).sign0;
      fz0 *= M_THIRD;//*(*ic).sign0;
      // Forces from k on j, factor M_THIRD for three-body
      fx1 *= M_THIRD;//*(*ic).sign1;
      fy1 *= M_THIRD;//*(*ic).sign1;
      fz1 *= M_THIRD;//*(*ic).sign1;
      // Figure out (i,j,k) indexes of particles present in ic
      long indexes[3], i;
      if ((*ic).sign0==1) {
        indexes[0] = (*delta0).j;
        indexes[1] = (*delta0).i;
      }
      else {
        indexes[1] = (*delta0).j;
        indexes[0] = (*delta0).i;
      }
      if ((*ic).sign1==1) {
        indexes[2] = (*delta1).j;
      }
      else {
        indexes[2] = (*delta1).i;
      }
      printf("INDEXES %5ld %5ld %5ld\n",indexes[0],indexes[1],indexes[2]);
      // Fill atomic virial tensors for every particle in ic
      for (i=0;i<3;i++) {
          avtens[indexes[i]*9  ] += fx0*(*delta0).dx+fx1*(*delta1).dx;
          avtens[indexes[i]*9+1] += fx0*(*delta0).dy+fx1*(*delta1).dy;
          avtens[indexes[i]*9+2] += fx0*(*delta0).dz+fx1*(*delta1).dz;
          avtens[indexes[i]*9+3] += fy0*(*delta0).dx+fy1*(*delta1).dx;
          avtens[indexes[i]*9+4] += fy0*(*delta0).dy+fy1*(*delta1).dy;
          avtens[indexes[i]*9+5] += fy0*(*delta0).dz+fy1*(*delta1).dz;
          avtens[indexes[i]*9+6] += fz0*(*delta0).dx+fz1*(*delta1).dx;
          avtens[indexes[i]*9+7] += fz0*(*delta0).dy+fz1*(*delta1).dy;
          avtens[indexes[i]*9+8] += fz0*(*delta0).dz+fz1*(*delta1).dz;
      }
  }
}

void back_bend_angle(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double *avtens) {
  back_bend_cos(ic, deltas, cos(value), -grad/sin(value), avtens);
}

void back_dihed_cos(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double *avtens) {
  long i;
  dlist_row_type *delta0, *delta1, *delta2;
  double a[3], b[3], dcos_da[3], dcos_db[3], da_ddel0[9], da_ddel1[9], db_ddel1[9];
  double dot0, dot2, n1, na, nb;
  double fx0, fx1, fx2, fy0, fy1, fy2, fz0, fz1, fz2;

  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  n1   = sqrt((*delta1).dx*(*delta1).dx + (*delta1).dy*(*delta1).dy + (*delta1).dz*(*delta1).dz);
  dot0 =      (*delta0).dx*(*delta1).dx + (*delta0).dy*(*delta1).dy + (*delta0).dz*(*delta1).dz;
  dot2 =      (*delta1).dx*(*delta2).dx + (*delta1).dy*(*delta2).dy + (*delta1).dz*(*delta2).dz;

  a[0] = ( (*delta0).dx - dot0*(*delta1).dx/(n1*n1) );
  a[1] = ( (*delta0).dy - dot0*(*delta1).dy/(n1*n1) );
  a[2] = ( (*delta0).dz - dot0*(*delta1).dz/(n1*n1) );
  b[0] = ( (*delta2).dx - dot2*(*delta1).dx/(n1*n1) );
  b[1] = ( (*delta2).dy - dot2*(*delta1).dy/(n1*n1) );
  b[2] = ( (*delta2).dz - dot2*(*delta1).dz/(n1*n1) );

  na = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  nb = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);

  value *= (*ic).sign0*(*ic).sign2;
  grad *= (*ic).sign0*(*ic).sign2;

  for (i=0; i<3; i++) {
    dcos_da[i] = (b[i]/nb - value*a[i]/na)/na;
    dcos_db[i] = (a[i]/na - value*b[i]/nb)/nb;
  }

  da_ddel0[0] = 1 - (*delta1).dx*(*delta1).dx/(n1*n1);
  da_ddel0[1] =   - (*delta1).dx*(*delta1).dy/(n1*n1);
  da_ddel0[2] =   - (*delta1).dx*(*delta1).dz/(n1*n1);
  da_ddel0[3] =   - (*delta1).dy*(*delta1).dx/(n1*n1);
  da_ddel0[4] = 1 - (*delta1).dy*(*delta1).dy/(n1*n1);
  da_ddel0[5] =   - (*delta1).dy*(*delta1).dz/(n1*n1);
  da_ddel0[6] =   - (*delta1).dz*(*delta1).dx/(n1*n1);
  da_ddel0[7] =   - (*delta1).dz*(*delta1).dy/(n1*n1);
  da_ddel0[8] = 1 - (*delta1).dz*(*delta1).dz/(n1*n1);

  da_ddel1[0] = ( - dot0/(n1*n1) - (*delta0).dx*(*delta1).dx/(n1*n1) + 2*dot0/(n1*n1*n1*n1)*(*delta1).dx*(*delta1).dx );
  da_ddel1[1] = (                - (*delta0).dx*(*delta1).dy/(n1*n1) + 2*dot0/(n1*n1*n1*n1)*(*delta1).dx*(*delta1).dy );
  da_ddel1[2] = (                - (*delta0).dx*(*delta1).dz/(n1*n1) + 2*dot0/(n1*n1*n1*n1)*(*delta1).dx*(*delta1).dz );
  da_ddel1[3] = (                - (*delta0).dy*(*delta1).dx/(n1*n1) + 2*dot0/(n1*n1*n1*n1)*(*delta1).dy*(*delta1).dx );
  da_ddel1[4] = ( - dot0/(n1*n1) - (*delta0).dy*(*delta1).dy/(n1*n1) + 2*dot0/(n1*n1*n1*n1)*(*delta1).dy*(*delta1).dy );
  da_ddel1[5] = (                - (*delta0).dy*(*delta1).dz/(n1*n1) + 2*dot0/(n1*n1*n1*n1)*(*delta1).dy*(*delta1).dz );
  da_ddel1[6] = (                - (*delta0).dz*(*delta1).dx/(n1*n1) + 2*dot0/(n1*n1*n1*n1)*(*delta1).dz*(*delta1).dx );
  da_ddel1[7] = (                - (*delta0).dz*(*delta1).dy/(n1*n1) + 2*dot0/(n1*n1*n1*n1)*(*delta1).dz*(*delta1).dy );
  da_ddel1[8] = ( - dot0/(n1*n1) - (*delta0).dz*(*delta1).dz/(n1*n1) + 2*dot0/(n1*n1*n1*n1)*(*delta1).dz*(*delta1).dz );

  db_ddel1[0] = ( - dot2/(n1*n1) - (*delta2).dx*(*delta1).dx/(n1*n1) + 2*dot2/(n1*n1*n1*n1)*(*delta1).dx*(*delta1).dx );
  db_ddel1[1] = (                - (*delta2).dx*(*delta1).dy/(n1*n1) + 2*dot2/(n1*n1*n1*n1)*(*delta1).dx*(*delta1).dy );
  db_ddel1[2] = (                - (*delta2).dx*(*delta1).dz/(n1*n1) + 2*dot2/(n1*n1*n1*n1)*(*delta1).dx*(*delta1).dz );
  db_ddel1[3] = (                - (*delta2).dy*(*delta1).dx/(n1*n1) + 2*dot2/(n1*n1*n1*n1)*(*delta1).dy*(*delta1).dx );
  db_ddel1[4] = ( - dot2/(n1*n1) - (*delta2).dy*(*delta1).dy/(n1*n1) + 2*dot2/(n1*n1*n1*n1)*(*delta1).dy*(*delta1).dy );
  db_ddel1[5] = (                - (*delta2).dy*(*delta1).dz/(n1*n1) + 2*dot2/(n1*n1*n1*n1)*(*delta1).dy*(*delta1).dz );
  db_ddel1[6] = (                - (*delta2).dz*(*delta1).dx/(n1*n1) + 2*dot2/(n1*n1*n1*n1)*(*delta1).dz*(*delta1).dx );
  db_ddel1[7] = (                - (*delta2).dz*(*delta1).dy/(n1*n1) + 2*dot2/(n1*n1*n1*n1)*(*delta1).dz*(*delta1).dy );
  db_ddel1[8] = ( - dot2/(n1*n1) - (*delta2).dz*(*delta1).dz/(n1*n1) + 2*dot2/(n1*n1*n1*n1)*(*delta1).dz*(*delta1).dz );

  fx0 = grad*(  dcos_da[0]*da_ddel0[0] + dcos_da[1]*da_ddel0[3] + dcos_da[2]*da_ddel0[6]);
  fy0 = grad*(  dcos_da[0]*da_ddel0[1] + dcos_da[1]*da_ddel0[4] + dcos_da[2]*da_ddel0[7]);
  fz0 = grad*(  dcos_da[0]*da_ddel0[2] + dcos_da[1]*da_ddel0[5] + dcos_da[2]*da_ddel0[8]);
  fx1 = grad*(  dcos_da[0]*da_ddel1[0] + dcos_da[1]*da_ddel1[3] + dcos_da[2]*da_ddel1[6]
                        + dcos_db[0]*db_ddel1[0] + dcos_db[1]*db_ddel1[3] + dcos_db[2]*db_ddel1[6]);
  fy1 = grad*(  dcos_da[0]*da_ddel1[1] + dcos_da[1]*da_ddel1[4] + dcos_da[2]*da_ddel1[7]
                        + dcos_db[0]*db_ddel1[1] + dcos_db[1]*db_ddel1[4] + dcos_db[2]*db_ddel1[7]);
  fz1 = grad*(  dcos_da[0]*da_ddel1[2] + dcos_da[1]*da_ddel1[5] + dcos_da[2]*da_ddel1[8]
                        + dcos_db[0]*db_ddel1[2] + dcos_db[1]*db_ddel1[5] + dcos_db[2]*db_ddel1[8]);
  fx2 = grad*(  dcos_db[0]*da_ddel0[0] + dcos_db[1]*da_ddel0[3] + dcos_db[2]*da_ddel0[6]);
  fy2 = grad*(  dcos_db[0]*da_ddel0[1] + dcos_db[1]*da_ddel0[4] + dcos_db[2]*da_ddel0[7]);
  fz2 = grad*(  dcos_db[0]*da_ddel0[2] + dcos_db[1]*da_ddel0[5] + dcos_db[2]*da_ddel0[8]);

  (*delta0).gx += fx0;
  (*delta0).gy += fy0;
  (*delta0).gz += fz0;
  (*delta1).gx += fx1;
  (*delta1).gy += fy1;
  (*delta1).gz += fz1;
  (*delta2).gx += fx2;
  (*delta2).gy += fy2;
  (*delta2).gz += fz2;

  if (avtens != NULL) {
      // Forces from i on j, factor 0.25 for four-body
      fx0 *= 0.25;
      fy0 *= 0.25;
      fz0 *= 0.25;
      // Forces from k on j, factor 0.25 for four-body
      fx1 *= 0.25;
      fy1 *= 0.25;
      fz1 *= 0.25;
      // Forces from l on k, factor 0.25 for four-body
      fx2 *= 0.25;
      fy2 *= 0.25;
      fz2 *= 0.25;
      // Figure out (i,j,k,l) indexes of particles present in ic
      long indexes[4], i;
      if ((*ic).sign0==1) {
        indexes[0] = (*delta0).j;
        indexes[1] = (*delta0).i;
      }
      else {
        indexes[1] = (*delta0).j;
        indexes[0] = (*delta0).i;
      }
      if ((*ic).sign1==1) {
        indexes[2] = (*delta1).j;
      }
      else {
        indexes[2] = (*delta1).i;
      }
      if ((*ic).sign2==1) {
        indexes[3] = (*delta2).j;
      }
      else {
        indexes[3] = (*delta2).i;
      }
      printf("INDEXES %5ld %5ld %5ld %5ld\n",indexes[0],indexes[1],indexes[2],indexes[3]);
      // Fill atomic virial tensors for every particle in ic
      for (i=0;i<4;i++) {
          avtens[indexes[i]*9  ] += fx0*(*delta0).dx+fx1*(*delta1).dx+fx2*(*delta2).dx;
          avtens[indexes[i]*9+1] += fx0*(*delta0).dy+fx1*(*delta1).dy+fx2*(*delta2).dy;
          avtens[indexes[i]*9+2] += fx0*(*delta0).dz+fx1*(*delta1).dz+fx2*(*delta2).dz;
          avtens[indexes[i]*9+3] += fy0*(*delta0).dx+fy1*(*delta1).dx+fy2*(*delta2).dx;
          avtens[indexes[i]*9+4] += fy0*(*delta0).dy+fy1*(*delta1).dy+fy2*(*delta2).dy;
          avtens[indexes[i]*9+5] += fy0*(*delta0).dz+fy1*(*delta1).dz+fy2*(*delta2).dz;
          avtens[indexes[i]*9+6] += fz0*(*delta0).dx+fz1*(*delta1).dx+fz2*(*delta2).dx;
          avtens[indexes[i]*9+7] += fz0*(*delta0).dy+fz1*(*delta1).dy+fz2*(*delta2).dy;
          avtens[indexes[i]*9+8] += fz0*(*delta0).dz+fz1*(*delta1).dz+fz2*(*delta2).dz;
      }
  }
}

void back_dihed_angle(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double *avtens) {
  double tmp = sin(value);
  if (tmp!=0.0) tmp = -grad/tmp;
  back_dihed_cos(ic, deltas, cos(value), tmp, avtens);
}

void back_oop_cos_low(iclist_row_type* ic, dlist_row_type *delta0, dlist_row_type *delta1, dlist_row_type *delta2, double value, double grad, double *avtens) {
  // This calculation is tedious. Expressions are checked with the following
  // maple commands (assuming the maple worksheet is bug-free)
  /*
  restart:
  with(LinearAlgebra):
  d0:=Vector([d0x,d0y,d0z]):
  d1:=Vector([d1x,d1y,d1z]):
  d2:=Vector([d2x,d2y,d2z]):
  assume(d0x,real,d0y,real,d0z,real,d1x,real,d1y,real,d1z,real,d2x,real,d2y,real,d2z,real):
  d0_cross_d1:=CrossProduct(d0,d1):
  d1_cross_d2:=CrossProduct(d1,d2):
  d2_cross_d0:=CrossProduct(d2,d0):
  n:=d0_cross_d1:
  n_sq:=Norm(n,2,conjugate=false)**2:
  d2_sq:=Norm(d2,2,conjugate=false)**2:
  n_dot_d2:=DotProduct(n,d2):
  fac:=n_dot_d2/d2_sq/n_sq:
  f:=DotProduct(d2,CrossProduct(d0,d1))/Norm(d2,2,conjugate=false)/Norm(CrossProduct(d0,d1),2,conjugate=false):
  cosphi:= (1-f**2)**(1/2):

  dcosphi_d0x:=diff(cosphi,d0x):
  dcosphi_d0y:=diff(cosphi,d0y):
  dcosphi_d0z:=diff(cosphi,d0z):
  dcosphi_d1x:=diff(cosphi,d1x):
  dcosphi_d1y:=diff(cosphi,d1y):
  dcosphi_d1z:=diff(cosphi,d1z):
  dcosphi_d2x:=diff(cosphi,d2x):
  dcosphi_d2y:=diff(cosphi,d2y):
  dcosphi_d2z:=diff(cosphi,d2z):

  simplify(dcosphi_d0x + fac/cosphi*( d1_cross_d2[1] - n_dot_d2/n_sq*(d1[2]*n[3]-d1[3]*n[2])));
  simplify(dcosphi_d0y + fac/cosphi*( d1_cross_d2[2] - n_dot_d2/n_sq*(d1[3]*n[1]-d1[1]*n[3])));
  simplify(dcosphi_d0z + fac/cosphi*( d1_cross_d2[3] - n_dot_d2/n_sq*(d1[1]*n[2]-d1[2]*n[1])));
  simplify(dcosphi_d1x + fac/cosphi*( d2_cross_d0[1] - n_dot_d2/n_sq*(d0[3]*n[2]-d0[2]*n[3])));
  simplify(dcosphi_d1y + fac/cosphi*( d2_cross_d0[2] - n_dot_d2/n_sq*(d0[1]*n[3]-d0[3]*n[1])));
  simplify(dcosphi_d1z + fac/cosphi*( d2_cross_d0[3] - n_dot_d2/n_sq*(d0[2]*n[1]-d0[1]*n[2])));
  simplify(dcosphi_d2x + fac/cosphi*( d0_cross_d1[1] - n_dot_d2/d2_sq*d2[1]) );
  simplify(dcosphi_d2y + fac/cosphi*( d0_cross_d1[2] - n_dot_d2/d2_sq*d2[2]) );
  simplify(dcosphi_d2z + fac/cosphi*( d0_cross_d1[3] - n_dot_d2/d2_sq*d2[3]) );
  */
  double n[3], d0_cross_d1[3], d1_cross_d2[3], d2_cross_d0[3];
  double n_sq, d2_sq, n_dot_d2, fac, tmp0, tmp1, tmp2;
  double fx0, fx1, fx2, fy0, fy1, fy2, fz0, fz1, fz2;
  // Cross products of delta vectors (introduce a function vectorproduct() ?)
  d0_cross_d1[0] = (*delta0).dy * (*delta1).dz - (*delta0).dz * (*delta1).dy;
  d0_cross_d1[1] = (*delta0).dz * (*delta1).dx - (*delta0).dx * (*delta1).dz;
  d0_cross_d1[2] = (*delta0).dx * (*delta1).dy - (*delta0).dy * (*delta1).dx;

  d1_cross_d2[0] = (*delta1).dy * (*delta2).dz - (*delta1).dz * (*delta2).dy;
  d1_cross_d2[1] = (*delta1).dz * (*delta2).dx - (*delta1).dx * (*delta2).dz;
  d1_cross_d2[2] = (*delta1).dx * (*delta2).dy - (*delta1).dy * (*delta2).dx;

  d2_cross_d0[0] = (*delta2).dy * (*delta0).dz - (*delta2).dz * (*delta0).dy;
  d2_cross_d0[1] = (*delta2).dz * (*delta0).dx - (*delta2).dx * (*delta0).dz;
  d2_cross_d0[2] = (*delta2).dx * (*delta0).dy - (*delta2).dy * (*delta0).dx;

  n[0] = d0_cross_d1[0], n[1] = d0_cross_d1[1],n[2] = d0_cross_d1[2]; //normal to plane of first two vectors
  // The squared norm of crossproduct of first two vectors
  n_sq = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  // Squared norm of the third vector
  d2_sq = (*delta2).dx * (*delta2).dx + (*delta2).dy * (*delta2).dy + (*delta2).dz * (*delta2).dz;
  // Dot product of the crossproduct of first two vectors with the third vector
  n_dot_d2 = n[0]*(*delta2).dx + n[1]*(*delta2).dy + n[2]*(*delta2).dz;
  // The expression for the cosine of the out-of-plane angle can be written as:
  // cos(phi) = sqrt(1-f**2), so the derivatives can be computed as
  // d cos(phi) / d x = -f / sqrt(1-f**2) * d f / d x
  fac = n_dot_d2/d2_sq/n_sq;
  tmp0 = fac/value;
  tmp1 = n_dot_d2/n_sq;
  tmp2 = n_dot_d2/d2_sq;
  fx0 = - tmp0*grad*( d1_cross_d2[0] -  tmp1*( (*delta1).dy*n[2] - (*delta1).dz*n[1] ) );
  fy0 = - tmp0*grad*( d1_cross_d2[1] -  tmp1*( (*delta1).dz*n[0] - (*delta1).dx*n[2] ) );
  fz0 = - tmp0*grad*( d1_cross_d2[2] -  tmp1*( (*delta1).dx*n[1] - (*delta1).dy*n[0] ) );
  fx1 = - tmp0*grad*( d2_cross_d0[0] -  tmp1*( (*delta0).dz*n[1] - (*delta0).dy*n[2] ) );
  fy1 = - tmp0*grad*( d2_cross_d0[1] -  tmp1*( (*delta0).dx*n[2] - (*delta0).dz*n[0] ) );
  fz1 = - tmp0*grad*( d2_cross_d0[2] -  tmp1*( (*delta0).dy*n[0] - (*delta0).dx*n[1] ) );
  fx2 = - tmp0*grad*( d0_cross_d1[0] -  tmp2*(*delta2).dx );
  fy2 = - tmp0*grad*( d0_cross_d1[1] -  tmp2*(*delta2).dy );
  fz2 = - tmp0*grad*( d0_cross_d1[2] -  tmp2*(*delta2).dz );

  (*delta0).gx += fx0;
  (*delta0).gy += fy0;
  (*delta0).gz += fz0;
  (*delta1).gx += fx1;
  (*delta1).gy += fy1;
  (*delta1).gz += fz1;
  (*delta2).gx += fx2;
  (*delta2).gy += fy2;
  (*delta2).gz += fz2;

  if (avtens != NULL) {
      // Forces from l on i, factor 0.25 for four-body
      fx0 *= 0.25;
      fy0 *= 0.25;
      fz0 *= 0.25;
      // Forces from l on j, factor 0.25 for four-body
      fx1 *= 0.25;
      fy1 *= 0.25;
      fz1 *= 0.25;
      // Forces from l on k, factor 0.25 for four-body
      fx2 *= 0.25;
      fy2 *= 0.25;
      fz2 *= 0.25;
      // Figure out (i,j,k,l) indexes of particles present in ic
      long indexes[4], i;
      if ((*ic).sign0==1) {
        indexes[0] = (*delta0).i;
        indexes[3] = (*delta0).j;
      }
      else {
        indexes[3] = (*delta0).i;
        indexes[0] = (*delta0).j;
      }
      if ((*ic).sign1==1) {
        indexes[1] = (*delta1).i;
      }
      else {
        indexes[1] = (*delta1).j;
      }
      if ((*ic).sign2==1) {
        indexes[2] = (*delta2).i;
      }
      else {
        indexes[2] = (*delta2).j;
      }
      printf("INDEXES %5ld %5ld %5ld %5ld\n",indexes[0],indexes[1],indexes[2],indexes[3]);
      printf("FORCES %+12.8f %+12.8f %+12.8f %+12.8f %+12.8f %+12.8f %+12.8f %+12.8f %+12.8f\n",fx0,fy0,fz0,fx1,fy1,fz1,fx2,fy2,fz2);
      printf("DISTS  %+12.8f %+12.8f %+12.8f %+12.8f %+12.8f %+12.8f %+12.8f %+12.8f %+12.8f\n",(*delta0).dx,(*delta0).dy,(*delta0).dz,(*delta1).dx,(*delta1).dy,(*delta1).dz,(*delta2).dx,(*delta2).dy,(*delta2).dz);
      // Fill atomic virial tensors for every particle in ic
      for (i=0;i<4;i++) {
          avtens[indexes[i]*9  ] += fx0*(*delta0).dx+fx1*(*delta1).dx+fx2*(*delta2).dx;
          avtens[indexes[i]*9+1] += fx0*(*delta0).dy+fx1*(*delta1).dy+fx2*(*delta2).dy;
          avtens[indexes[i]*9+2] += fx0*(*delta0).dz+fx1*(*delta1).dz+fx2*(*delta2).dz;
          avtens[indexes[i]*9+3] += fy0*(*delta0).dx+fy1*(*delta1).dx+fy2*(*delta2).dx;
          avtens[indexes[i]*9+4] += fy0*(*delta0).dy+fy1*(*delta1).dy+fy2*(*delta2).dy;
          avtens[indexes[i]*9+5] += fy0*(*delta0).dz+fy1*(*delta1).dz+fy2*(*delta2).dz;
          avtens[indexes[i]*9+6] += fz0*(*delta0).dx+fz1*(*delta1).dx+fz2*(*delta2).dx;
          avtens[indexes[i]*9+7] += fz0*(*delta0).dy+fz1*(*delta1).dy+fz2*(*delta2).dy;
          avtens[indexes[i]*9+8] += fz0*(*delta0).dz+fz1*(*delta1).dz+fz2*(*delta2).dz;
      }
  }
}

void back_oop_cos(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double *avtens) {
  dlist_row_type *delta0, *delta1, *delta2;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  back_oop_cos_low(ic, delta0, delta1, delta2, value, grad, avtens);
}

void back_oop_meancos(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double *avtens) {
  dlist_row_type *delta0, *delta1, *delta2;
  double tmp;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  tmp = forward_oop_cos_low((double*)delta0, (double*)delta1, (double*)delta2);
  back_oop_cos_low(ic, delta0, delta1, delta2, tmp, grad/3.0, avtens);
  tmp = forward_oop_cos_low((double*)delta2, (double*)delta0,(double*)delta1);
  back_oop_cos_low(ic, delta2, delta0, delta1, tmp, grad/3.0, avtens);
  tmp = forward_oop_cos_low((double*)delta1, (double*)delta2, (double*)delta0);
  back_oop_cos_low(ic, delta1, delta2, delta0, tmp, grad/3.0, avtens);
}

void back_oop_angle_low(iclist_row_type* ic, dlist_row_type *delta0, dlist_row_type *delta1, dlist_row_type *delta2, double value, double grad, double *avtens) {
  double tmp = sin(value);
  if (tmp!=0.0) tmp = -grad/tmp;
  back_oop_cos_low(ic, delta0, delta1, delta2, cos(value), tmp, avtens);
}

void back_oop_angle(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double *avtens) {
  dlist_row_type *delta0, *delta1, *delta2;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  back_oop_angle_low(ic, delta0, delta1, delta2, value, grad, avtens);
}

void back_oop_meanangle(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double *avtens) {
  dlist_row_type *delta0, *delta1, *delta2;
  double tmp;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  tmp = forward_oop_angle_low((double*)delta0, (double*)delta1, (double*)delta2);
  back_oop_angle_low(ic, delta0, delta1, delta2, tmp, grad/3.0, avtens);
  tmp = forward_oop_angle_low((double*)delta2, (double*)delta0, (double*)delta1);
  back_oop_angle_low(ic, delta2, delta0, delta1, tmp, grad/3.0, avtens);
  tmp = forward_oop_angle_low((double*)delta1, (double*)delta2, (double*)delta0);
  back_oop_angle_low(ic, delta1, delta2, delta0, tmp, grad/3.0, avtens);
}

void back_oop_distance(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double *avtens) {
  dlist_row_type *delta0, *delta1, *delta2;
  double n[3], d0_cross_d1[3], d1_cross_d2[3], d2_cross_d0[3];
  double n_norm, n_dot_d2, fac, tmp0;
  double fx0, fx1, fx2, fy0, fy1, fy2, fz0, fz1, fz2;

  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  // Cross products of delta vectors (introduce a function vectorproduct() ?)
  d0_cross_d1[0] = (*delta0).dy * (*delta1).dz - (*delta0).dz * (*delta1).dy;
  d0_cross_d1[1] = (*delta0).dz * (*delta1).dx - (*delta0).dx * (*delta1).dz;
  d0_cross_d1[2] = (*delta0).dx * (*delta1).dy - (*delta0).dy * (*delta1).dx;

  d1_cross_d2[0] = (*delta1).dy * (*delta2).dz - (*delta1).dz * (*delta2).dy;
  d1_cross_d2[1] = (*delta1).dz * (*delta2).dx - (*delta1).dx * (*delta2).dz;
  d1_cross_d2[2] = (*delta1).dx * (*delta2).dy - (*delta1).dy * (*delta2).dx;

  d2_cross_d0[0] = (*delta2).dy * (*delta0).dz - (*delta2).dz * (*delta0).dy;
  d2_cross_d0[1] = (*delta2).dz * (*delta0).dx - (*delta2).dx * (*delta0).dz;
  d2_cross_d0[2] = (*delta2).dx * (*delta0).dy - (*delta2).dy * (*delta0).dx;

  n[0] = d0_cross_d1[0], n[1] = d0_cross_d1[1],n[2] = d0_cross_d1[2]; //normal to plane of first two vectors
  // The squared norm of crossproduct of first two vectors
  n_norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  // If n_norm==0, the first and second vector don't span a plane!
  //if (n_norm == 0) ???;
  // Dot product of the crossproduct of first two vectors with the third vector
  n_dot_d2 = n[0]*(*delta2).dx + n[1]*(*delta2).dy + n[2]*(*delta2).dz;
  fac = grad/n_norm*(*ic).sign0*(*ic).sign1*(*ic).sign2;
  tmp0 = n_dot_d2/n_norm/n_norm;

  fx0 = fac*( d1_cross_d2[0] -  tmp0*( (*delta1).dy*n[2] - (*delta1).dz*n[1] ) );
  fy0 = fac*( d1_cross_d2[1] -  tmp0*( (*delta1).dz*n[0] - (*delta1).dx*n[2] ) );
  fz0 = fac*( d1_cross_d2[2] -  tmp0*( (*delta1).dx*n[1] - (*delta1).dy*n[0] ) );
  fx1 = fac*( d2_cross_d0[0] -  tmp0*( (*delta0).dz*n[1] - (*delta0).dy*n[2] ) );
  fy1 = fac*( d2_cross_d0[1] -  tmp0*( (*delta0).dx*n[2] - (*delta0).dz*n[0] ) );
  fz1 = fac*( d2_cross_d0[2] -  tmp0*( (*delta0).dy*n[0] - (*delta0).dx*n[1] ) );
  fx2 = fac*( d0_cross_d1[0] );
  fy2 = fac*( d0_cross_d1[1] );
  fz2 = fac*( d0_cross_d1[2] );

  (*delta0).gx += fx0;
  (*delta0).gy += fy0;
  (*delta0).gz += fz0;
  (*delta1).gx += fx1;
  (*delta1).gy += fy1;
  (*delta1).gz += fz1;
  (*delta2).gx += fx2;
  (*delta2).gy += fy2;
  (*delta2).gz += fz2;

  if (avtens != NULL) {
      // Forces from l on i, factor 0.25 for four-body
      fx0 *= 0.25;
      fy0 *= 0.25;
      fz0 *= 0.25;
      // Forces from l on j, factor 0.25 for four-body
      fx1 *= 0.25;
      fy1 *= 0.25;
      fz1 *= 0.25;
      // Forces from l on k, factor 0.25 for four-body
      fx2 *= 0.25;
      fy2 *= 0.25;
      fz2 *= 0.25;
      // Figure out (i,j,k,l) indexes of particles present in ic
      long indexes[4], i;
      if ((*ic).sign0==1) {
        indexes[0] = (*delta0).i;
        indexes[1] = (*delta0).j;
      }
      else {
        indexes[1] = (*delta0).i;
        indexes[0] = (*delta0).j;
      }
      if ((*ic).sign1==1) {
        indexes[2] = (*delta1).i;
      }
      else {
        indexes[2] = (*delta1).j;
      }
      if ((*ic).sign2==1) {
        indexes[3] = (*delta2).i;
      }
      else {
        indexes[3] = (*delta2).j;
      }
      printf("INDEXES %5ld %5ld %5ld %5ld\n",indexes[0],indexes[1],indexes[2],indexes[3]);
      // Fill atomic virial tensors for every particle in ic
      for (i=0;i<4;i++) {
          avtens[indexes[i]*9  ] += fx0*(*delta0).dx+fx1*(*delta1).dx+fx2*(*delta2).dx;
          avtens[indexes[i]*9+1] += fx0*(*delta0).dy+fx1*(*delta1).dy+fx2*(*delta2).dy;
          avtens[indexes[i]*9+2] += fx0*(*delta0).dz+fx1*(*delta1).dz+fx2*(*delta2).dz;
          avtens[indexes[i]*9+3] += fy0*(*delta0).dx+fy1*(*delta1).dx+fy2*(*delta2).dx;
          avtens[indexes[i]*9+4] += fy0*(*delta0).dy+fy1*(*delta1).dy+fy2*(*delta2).dy;
          avtens[indexes[i]*9+5] += fy0*(*delta0).dz+fy1*(*delta1).dz+fy2*(*delta2).dz;
          avtens[indexes[i]*9+6] += fz0*(*delta0).dx+fz1*(*delta1).dx+fz2*(*delta2).dx;
          avtens[indexes[i]*9+7] += fz0*(*delta0).dy+fz1*(*delta1).dy+fz2*(*delta2).dy;
          avtens[indexes[i]*9+8] += fz0*(*delta0).dz+fz1*(*delta1).dz+fz2*(*delta2).dz;
      }
  }
}

ic_back_type ic_back_fns[11] = {
  back_bond, back_bend_cos, back_bend_angle, back_dihed_cos, back_dihed_angle, back_bond,
  back_oop_cos, back_oop_meancos, back_oop_angle, back_oop_meanangle, back_oop_distance
};

void iclist_back(dlist_row_type* deltas, iclist_row_type* ictab, long nic, double *avtens) {
  long i;
  for (i=0; i<nic; i++) {
    ic_back_fns[ictab[i].kind](ictab + i, deltas, ictab[i].value, ictab[i].grad, avtens);
  }
}
