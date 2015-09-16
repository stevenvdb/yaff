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
  //TODO Check what we should do with (*ic).signs here
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
  return n_dot_d2/n_norm;
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


typedef void (*ic_back_type)(iclist_row_type*, dlist_row_type*, double, double);

void back_bond(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad) {
  dlist_row_type *delta;
  double x;
  delta = deltas + (*ic).i0;
  x = grad/value;
  (*delta).gx += x*(*delta).dx;
  (*delta).gy += x*(*delta).dy;
  (*delta).gz += x*(*delta).dz;
}

void back_bend_cos(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad) {
  dlist_row_type *delta0, *delta1;
  double e0[3], e1[3];
  double d0, d1, fac;
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
  (*delta0).gx += fac*(e1[0] - value*e0[0]);
  (*delta0).gy += fac*(e1[1] - value*e0[1]);
  (*delta0).gz += fac*(e1[2] - value*e0[2]);
  fac = grad/d1;
  (*delta1).gx += fac*(e0[0] - value*e1[0]);
  (*delta1).gy += fac*(e0[1] - value*e1[1]);
  (*delta1).gz += fac*(e0[2] - value*e1[2]);
}

void back_bend_angle(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad) {
  back_bend_cos(ic, deltas, cos(value), -grad/sin(value));
}

void back_dihed_cos(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad) {
  long i;
  dlist_row_type *delta0, *delta1, *delta2;
  double a[3], b[3], dcos_da[3], dcos_db[3], da_ddel0[9], da_ddel1[9], db_ddel1[9];
  double dot0, dot2, n1, na, nb;

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

  (*delta0).gx += grad*(  dcos_da[0]*da_ddel0[0] + dcos_da[1]*da_ddel0[3] + dcos_da[2]*da_ddel0[6]);
  (*delta0).gy += grad*(  dcos_da[0]*da_ddel0[1] + dcos_da[1]*da_ddel0[4] + dcos_da[2]*da_ddel0[7]);
  (*delta0).gz += grad*(  dcos_da[0]*da_ddel0[2] + dcos_da[1]*da_ddel0[5] + dcos_da[2]*da_ddel0[8]);
  (*delta1).gx += grad*(  dcos_da[0]*da_ddel1[0] + dcos_da[1]*da_ddel1[3] + dcos_da[2]*da_ddel1[6]
                        + dcos_db[0]*db_ddel1[0] + dcos_db[1]*db_ddel1[3] + dcos_db[2]*db_ddel1[6]);
  (*delta1).gy += grad*(  dcos_da[0]*da_ddel1[1] + dcos_da[1]*da_ddel1[4] + dcos_da[2]*da_ddel1[7]
                        + dcos_db[0]*db_ddel1[1] + dcos_db[1]*db_ddel1[4] + dcos_db[2]*db_ddel1[7]);
  (*delta1).gz += grad*(  dcos_da[0]*da_ddel1[2] + dcos_da[1]*da_ddel1[5] + dcos_da[2]*da_ddel1[8]
                        + dcos_db[0]*db_ddel1[2] + dcos_db[1]*db_ddel1[5] + dcos_db[2]*db_ddel1[8]);
  (*delta2).gx += grad*(  dcos_db[0]*da_ddel0[0] + dcos_db[1]*da_ddel0[3] + dcos_db[2]*da_ddel0[6]);
  (*delta2).gy += grad*(  dcos_db[0]*da_ddel0[1] + dcos_db[1]*da_ddel0[4] + dcos_db[2]*da_ddel0[7]);
  (*delta2).gz += grad*(  dcos_db[0]*da_ddel0[2] + dcos_db[1]*da_ddel0[5] + dcos_db[2]*da_ddel0[8]);
}

void back_dihed_angle(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad) {
  double tmp = sin(value);
  if (tmp!=0.0) tmp = -grad/tmp;
  back_dihed_cos(ic, deltas, cos(value), tmp);
}

void back_oop_cos_low(dlist_row_type *delta0, dlist_row_type *delta1, dlist_row_type *delta2, double value, double grad) {
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
  (*delta0).gx += - tmp0*grad*( d1_cross_d2[0] -  tmp1*( (*delta1).dy*n[2] - (*delta1).dz*n[1] ) );
  (*delta0).gy += - tmp0*grad*( d1_cross_d2[1] -  tmp1*( (*delta1).dz*n[0] - (*delta1).dx*n[2] ) );
  (*delta0).gz += - tmp0*grad*( d1_cross_d2[2] -  tmp1*( (*delta1).dx*n[1] - (*delta1).dy*n[0] ) );
  (*delta1).gx += - tmp0*grad*( d2_cross_d0[0] -  tmp1*( (*delta0).dz*n[1] - (*delta0).dy*n[2] ) );
  (*delta1).gy += - tmp0*grad*( d2_cross_d0[1] -  tmp1*( (*delta0).dx*n[2] - (*delta0).dz*n[0] ) );
  (*delta1).gz += - tmp0*grad*( d2_cross_d0[2] -  tmp1*( (*delta0).dy*n[0] - (*delta0).dx*n[1] ) );
  (*delta2).gx += - tmp0*grad*( d0_cross_d1[0] -  tmp2*(*delta2).dx );
  (*delta2).gy += - tmp0*grad*( d0_cross_d1[1] -  tmp2*(*delta2).dy );
  (*delta2).gz += - tmp0*grad*( d0_cross_d1[2] -  tmp2*(*delta2).dz );
}

void back_oop_cos(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad) {
  dlist_row_type *delta0, *delta1, *delta2;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  back_oop_cos_low(delta0, delta1, delta2, value, grad);
}

void back_oop_meancos(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad) {
  dlist_row_type *delta0, *delta1, *delta2;
  double tmp;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  tmp = forward_oop_cos_low((double*)delta0, (double*)delta1, (double*)delta2);
  back_oop_cos_low(delta0, delta1, delta2, tmp, grad/3.0);
  tmp = forward_oop_cos_low((double*)delta2, (double*)delta0,(double*)delta1);
  back_oop_cos_low(delta2, delta0, delta1, tmp, grad/3.0);
  tmp = forward_oop_cos_low((double*)delta1, (double*)delta2, (double*)delta0);
  back_oop_cos_low(delta1, delta2, delta0, tmp, grad/3.0);
}

void back_oop_angle_low(dlist_row_type *delta0, dlist_row_type *delta1, dlist_row_type *delta2, double value, double grad) {
  double tmp = sin(value);
  if (tmp!=0.0) tmp = -grad/tmp;
  back_oop_cos_low(delta0, delta1, delta2, cos(value), tmp);
}

void back_oop_angle(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad) {
  dlist_row_type *delta0, *delta1, *delta2;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  back_oop_angle_low(delta0, delta1, delta2, value, grad);
}

void back_oop_meanangle(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad) {
  dlist_row_type *delta0, *delta1, *delta2;
  double tmp;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  tmp = forward_oop_angle_low((double*)delta0, (double*)delta1, (double*)delta2);
  back_oop_angle_low(delta0, delta1, delta2, tmp, grad/3.0);
  tmp = forward_oop_angle_low((double*)delta2, (double*)delta0, (double*)delta1);
  back_oop_angle_low(delta2, delta0, delta1, tmp, grad/3.0);
  tmp = forward_oop_angle_low((double*)delta1, (double*)delta2, (double*)delta0);
  back_oop_angle_low(delta1, delta2, delta0, tmp, grad/3.0);
}

void back_oop_distance(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad) {
  dlist_row_type *delta0, *delta1, *delta2;
  double n[3], d0_cross_d1[3], d1_cross_d2[3], d2_cross_d0[3];
  double n_norm, n_dot_d2, fac, tmp0;

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
  fac = grad/n_norm;
  tmp0 = n_dot_d2/n_norm/n_norm;

  (*delta0).gx += fac*( d1_cross_d2[0] -  tmp0*( (*delta1).dy*n[2] - (*delta1).dz*n[1] ) );
  (*delta0).gy += fac*( d1_cross_d2[1] -  tmp0*( (*delta1).dz*n[0] - (*delta1).dx*n[2] ) );
  (*delta0).gz += fac*( d1_cross_d2[2] -  tmp0*( (*delta1).dx*n[1] - (*delta1).dy*n[0] ) );
  (*delta1).gx += fac*( d2_cross_d0[0] -  tmp0*( (*delta0).dz*n[1] - (*delta0).dy*n[2] ) );
  (*delta1).gy += fac*( d2_cross_d0[1] -  tmp0*( (*delta0).dx*n[2] - (*delta0).dz*n[0] ) );
  (*delta1).gz += fac*( d2_cross_d0[2] -  tmp0*( (*delta0).dy*n[0] - (*delta0).dx*n[1] ) );
  (*delta2).gx += fac*( d0_cross_d1[0] );
  (*delta2).gy += fac*( d0_cross_d1[1] );
  (*delta2).gz += fac*( d0_cross_d1[2] );
}

ic_back_type ic_back_fns[11] = {
  back_bond, back_bend_cos, back_bend_angle, back_dihed_cos, back_dihed_angle, back_bond,
  back_oop_cos, back_oop_meancos, back_oop_angle, back_oop_meanangle, back_oop_distance
};

void iclist_back(dlist_row_type* deltas, iclist_row_type* ictab, long nic) {
  long i;
  for (i=0; i<nic; i++) {
    ic_back_fns[ictab[i].kind](ictab + i, deltas, ictab[i].value, ictab[i].grad);
  }
}

typedef void (*ic_jacobian_type)(iclist_row_type*, dlist_row_type*, double, double *, long, long);

void jacobian_1delta(iclist_row_type* ic, dlist_row_type* deltas, double value, double* jacobian, long iic, long ndelta) {
  dlist_row_type *delta;
  delta = deltas + (*ic).i0;
  (*delta).gx = 0.0;
  (*delta).gy = 0.0;
  (*delta).gz = 0.0;
  ic_back_fns[(*ic).kind](ic, deltas, value, 1.0);
  jacobian[iic*ndelta+3*(*ic).i0+0] += (*delta).gx;
  jacobian[iic*ndelta+3*(*ic).i0+1] += (*delta).gy;
  jacobian[iic*ndelta+3*(*ic).i0+2] += (*delta).gz;
}

void jacobian_2delta(iclist_row_type* ic, dlist_row_type* deltas, double value, double* jacobian, long iic, long ndelta) {
  dlist_row_type *delta0, *delta1;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  (*delta0).gx = 0.0;
  (*delta0).gy = 0.0;
  (*delta0).gz = 0.0;
  (*delta1).gx = 0.0;
  (*delta1).gy = 0.0;
  (*delta1).gz = 0.0;
  ic_back_fns[(*ic).kind](ic, deltas, value, 1.0);
  jacobian[iic*ndelta+3*(*ic).i0+0] += (*delta0).gx;
  jacobian[iic*ndelta+3*(*ic).i0+1] += (*delta0).gy;
  jacobian[iic*ndelta+3*(*ic).i0+2] += (*delta0).gz;
  jacobian[iic*ndelta+3*(*ic).i1+0] += (*delta1).gx;
  jacobian[iic*ndelta+3*(*ic).i1+1] += (*delta1).gy;
  jacobian[iic*ndelta+3*(*ic).i1+2] += (*delta1).gz;
}

void jacobian_3delta(iclist_row_type* ic, dlist_row_type* deltas, double value, double* jacobian, long iic, long ndelta) {
  dlist_row_type *delta0, *delta1, *delta2;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  (*delta0).gx = 0.0;
  (*delta0).gy = 0.0;
  (*delta0).gz = 0.0;
  (*delta1).gx = 0.0;
  (*delta1).gy = 0.0;
  (*delta1).gz = 0.0;
  (*delta2).gx = 0.0;
  (*delta2).gy = 0.0;
  (*delta2).gz = 0.0;
  ic_back_fns[(*ic).kind](ic, deltas, value, 1.0);
  jacobian[iic*ndelta+3*(*ic).i0+0] += (*delta0).gx;
  jacobian[iic*ndelta+3*(*ic).i0+1] += (*delta0).gy;
  jacobian[iic*ndelta+3*(*ic).i0+2] += (*delta0).gz;
  jacobian[iic*ndelta+3*(*ic).i1+0] += (*delta1).gx;
  jacobian[iic*ndelta+3*(*ic).i1+1] += (*delta1).gy;
  jacobian[iic*ndelta+3*(*ic).i1+2] += (*delta1).gz;
  jacobian[iic*ndelta+3*(*ic).i2+0] += (*delta2).gx;
  jacobian[iic*ndelta+3*(*ic).i2+1] += (*delta2).gy;
  jacobian[iic*ndelta+3*(*ic).i2+2] += (*delta2).gz;
}

ic_jacobian_type ic_jacobian_fns[11] = {
  jacobian_1delta, jacobian_2delta, jacobian_2delta, jacobian_3delta, jacobian_3delta, jacobian_1delta,
  jacobian_3delta, jacobian_3delta, jacobian_3delta, jacobian_3delta, jacobian_3delta
};

void iclist_jacobian(dlist_row_type* deltas, iclist_row_type* ictab, long nic, long ndelta, double* jacobian) {
  long i;
  for (i=0; i<nic; i++) {
    ic_jacobian_fns[ictab[i].kind](ictab + i, deltas, ictab[i].value, jacobian, i, ndelta);
  }
}

typedef void (*ic_hessian_type)(iclist_row_type*, dlist_row_type*, double, double, double *, long);

void hessian_bond(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double* hessian, long ndelta) {
   // Initialization
   dlist_row_type *delta0;
   delta0 = deltas + (*ic).i0;
   double x0 = (*delta0).dx;
   double y0 = (*delta0).dy;
   double z0 = (*delta0).dz;
   // Compute auxiliary variables
   double t0 = y0*y0;
   double t1 = z0*z0;
   double t2 = x0*x0;
   double t3 = t0 + t2;
   double t4 = 1.0/(sqrt((t1 + t3)*(t1 + t3)*(t1 + t3)));
   double t5 = grad*t4;
   double t6 = grad*t4*x0;
   double t7 = -t6*y0;
   double t8 = -t6*z0;
   double t9 = -t5*y0*z0;
   // Fill up Hessian
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i0+0] =  t5*(t0 + t1);
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i0+1] =  t7;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i0+2] =  t8;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i0+0] =  t7;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i0+1] =  t5*(t1 + t2);
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i0+2] =  t9;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i0+0] =  t8;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i0+1] =  t9;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i0+2] =  t3*t5;
}


void hessian_bend_cos(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double* hessian, long ndelta) {
   // Initialization
   dlist_row_type *delta0;
   delta0 = deltas + (*ic).i0;
   double x0 = (*delta0).dx;
   double y0 = (*delta0).dy;
   double z0 = (*delta0).dz;
   dlist_row_type *delta1;
   delta1 = deltas + (*ic).i1;
   double x1 = (*delta1).dx;
   double y1 = (*delta1).dy;
   double z1 = (*delta1).dz;
   double s0 = (*ic).sign0;
   double s1 = (*ic).sign1;
   // Compute auxiliary variables
   double t0 = x0*x0;
   double t1 = y0*y0;
   double t2 = z0*z0;
   double t3 = t0 + t1 + t2;
   double t4 = x1*x1;
   double t5 = y1*y1;
   double t6 = z1*z1;
   double t7 = t4 + t5 + t6;
   double t8 = grad*s0*s1/(sqrt(t7)*sqrt(pow(t3, 5)));
   double t9 = x0*x1;
   double t10 = 3*t9;
   double t11 = y0*y1;
   double t12 = 3*t11;
   double t13 = z0*z1;
   double t14 = 3*t13;
   double t15 = t10 + t12 + t14;
   double t16 = t10 + t11 + t13;
   double t17 = t11 + t9;
   double t18 = t13 + t17;
   double t19 = 3*t18*x0;
   double t20 = x0*y1;
   double t21 = x1*y0;
   double t22 = t20 + t21;
   double t23 = t8*(t19*y0 - t22*t3);
   double t24 = x0*z1;
   double t25 = x1*z0;
   double t26 = t24 + t25;
   double t27 = t8*(t19*z0 - t26*t3);
   double t28 = grad*s0*s1/(sqrt(pow(t3, 3))*sqrt(pow(t7, 3)));
   double t29 = t3*t7;
   double t30 = t28*(-t0*t7 + t18*t9 + t29 - t3*t4);
   double t31 = t7*x0;
   double t32 = t3*x1;
   double t33 = -t31*y0 - t32*y1;
   double t34 = t28*(t18*t20 + t33);
   double t35 = -t31*z0 - t32*z1;
   double t36 = t28*(t18*t24 + t35);
   double t37 = t12 + t13 + t9;
   double t38 = y0*z0;
   double t39 = y0*z1;
   double t40 = y1*z0;
   double t41 = t39 + t40;
   double t42 = t8*(t15*t38 - t3*t41);
   double t43 = t28*(t18*t21 + t33);
   double t44 = t28*(-t1*t7 + t11*t18 + t29 - t3*t5);
   double t45 = y1*z1;
   double t46 = -t3*t45 - t38*t7;
   double t47 = t28*(t18*t39 + t46);
   double t48 = t14 + t17;
   double t49 = t28*(t18*t25 + t35);
   double t50 = t28*(t18*t40 + t46);
   double t51 = t28*(t13*t18 - t2*t7 + t29 - t3*t6);
   double t52 = grad*s0*s1/(sqrt(t3)*sqrt(pow(t7, 5)));
   double t53 = 3*t18*x1;
   double t54 = t52*(-t22*t7 + t53*y1);
   double t55 = t52*(-t26*t7 + t53*z1);
   double t56 = t52*(t15*t45 - t41*t7);
   // Fill up Hessian
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i0+0] +=  t8*(t0*t15 - t16*t3);
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i0+1] +=  t23;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i0+2] +=  t27;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i0+0] +=  t23;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i0+1] +=  t8*(t1*t15 - t3*t37);
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i0+2] +=  t42;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i0+0] +=  t27;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i0+1] +=  t42;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i0+2] +=  t8*(t15*t2 - t3*t48);
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i1+0] +=  t30;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i1+1] +=  t34;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i1+2] +=  t36;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i1+0] +=  t43;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i1+1] +=  t44;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i1+2] +=  t47;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i1+0] +=  t49;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i1+1] +=  t50;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i1+2] +=  t51;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i0+0] +=  t30;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i0+1] +=  t43;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i0+2] +=  t49;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i0+0] +=  t34;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i0+1] +=  t44;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i0+2] +=  t50;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i0+0] +=  t36;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i0+1] +=  t47;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i0+2] +=  t51;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i1+0] +=  t52*(t15*t4 - t16*t7);
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i1+1] +=  t54;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i1+2] +=  t55;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i1+0] +=  t54;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i1+1] +=  t52*(t15*t5 - t37*t7);
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i1+2] +=  t56;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i1+0] +=  t55;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i1+1] +=  t56;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i1+2] +=  t52*(t15*t6 - t48*t7);
}


void hessian_bend_angle(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double* hessian, long ndelta) {
   dlist_row_type *delta0, *delta1;
   // First part
   hessian_bend_cos(ic, deltas, cos(value), -grad/sin(value), hessian, ndelta);
   // Second part
   double jacobian[6];
   delta0 = deltas + (*ic).i0;
   delta1 = deltas + (*ic).i1;
   (*delta0).gx = 0.0;
   (*delta0).gy = 0.0;
   (*delta0).gz = 0.0;
   (*delta1).gx = 0.0;
   (*delta1).gy = 0.0;
   (*delta1).gz = 0.0;
   ic_back_fns[(*ic).kind](ic, deltas, value, 1.0);
   jacobian[0] = (*delta0).gx;
   jacobian[1] = (*delta0).gy;
   jacobian[2] = (*delta0).gz;
   jacobian[3] = (*delta1).gx;
   jacobian[4] = (*delta1).gy;
   jacobian[5] = (*delta1).gz;
   long i,j;
   double fac;
   fac = grad*cos(value)/sin(value);
   for (i=0;i<3;i++){
     for (j=0;j<3;j++){
       hessian[(3*(*ic).i0+i)*ndelta + 3*(*ic).i0+j] -= fac*jacobian[i]*jacobian[j];
       hessian[(3*(*ic).i0+i)*ndelta + 3*(*ic).i1+j] -= fac*jacobian[i]*jacobian[j+3];
       hessian[(3*(*ic).i1+i)*ndelta + 3*(*ic).i0+j] -= fac*jacobian[i+3]*jacobian[j];
       hessian[(3*(*ic).i1+i)*ndelta + 3*(*ic).i1+j] -= fac*jacobian[i+3]*jacobian[j+3];
     }
   }
}

void hessian_dihed_cos(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double* hessian, long ndelta) {
   // Initialization
   dlist_row_type *delta0;
   delta0 = deltas + (*ic).i0;
   double x0 = (*delta0).dx;
   double y0 = (*delta0).dy;
   double z0 = (*delta0).dz;
   dlist_row_type *delta1;
   delta1 = deltas + (*ic).i1;
   double x1 = (*delta1).dx;
   double y1 = (*delta1).dy;
   double z1 = (*delta1).dz;
   dlist_row_type *delta2;
   delta2 = deltas + (*ic).i2;
   double x2 = (*delta2).dx;
   double y2 = (*delta2).dy;
   double z2 = (*delta2).dz;
   double s0 = (*ic).sign0;
   double s2 = (*ic).sign2;
   // Compute auxiliary variables
   double t0 = x1*x1;
   double t1 = y1*y1;
   double t2 = z1*z1;
   double t3 = t0 + t1 + t2;
   double t4 = 1.0/t3;
   double t5 = t0*t4;
   double t6 = -t5 + 1;
   double t7 = 1.0/(t3*t3);
   double t8 = t0*t7;
   double t9 = t1*t8;
   double t10 = -t9;
   double t11 = t2*t8;
   double t12 = -t11;
   double t13 = t10 + t12;
   double t14 = x0*x1;
   double t15 = y0*y1;
   double t16 = z0*z1;
   double t17 = t14 + t15 + t16;
   double t18 = t4*x1;
   double t19 = -t17*t18 + x0;
   double t20 = t4*y1;
   double t21 = -t17*t20 + y0;
   double t22 = t4*z1;
   double t23 = -t17*t22 + z0;
   double t24 = t19*t19 + t21*t21 + t23*t23;
   double t25 = 1.0/(sqrt(t24*t24*t24));
   double t26 = x1*x2;
   double t27 = y1*y2;
   double t28 = z1*z2;
   double t29 = t26 + t27 + t28;
   double t30 = t18*t29;
   double t31 = t30 - x2;
   double t32 = t20*t29;
   double t33 = t32 - y2;
   double t34 = t22*t29;
   double t35 = t34 - z2;
   double t36 = t31*t31 + t33*t33 + t35*t35;
   double t37 = 1.0/(sqrt(t36));
   double t38 = -t30 + x2;
   double t39 = -t32 + y2;
   double t40 = -t34 + z2;
   double t41 = t19*t38 + t21*t39 + t23*t40;
   double t42 = grad*s0*s2*t25*t37*t41;
   double t43 = t4*x1*y1;
   double t44 = t21*t43;
   double t45 = t4*x1*z1;
   double t46 = t23*t45;
   double t47 = -t19*t6 + t44 + t46;
   double t48 = grad*s0*s2*t25*t37*t47;
   double t49 = 2*t5;
   double t50 = -t49 + 2;
   double t51 = -3.0L/2.0L*t17*t4*x1 + (3.0L/2.0L)*x0;
   double t52 = 1.0/(sqrt(t24*t24*t24*t24*t24));
   double t53 = grad*s0*s2*t37*t41*t47*t52;
   double t54 = t2*t7*x1;
   double t55 = t54*y1;
   double t56 = -t55;
   double t57 = t1*t4;
   double t58 = -t57 + 1;
   double t59 = t43*t58;
   double t60 = (1.0L/2.0L)*t4*x1*y1;
   double t61 = t38*t6 - t39*t43 - t40*t45;
   double t62 = t19*t43;
   double t63 = t4*y1*z1;
   double t64 = t23*t63;
   double t65 = -t21*t58 + t62 + t64;
   double t66 = grad*s0*s2*t25*t37*t65;
   double t67 = -t38*t43 + t39*t58 - t40*t63;
   double t68 = 2*t57;
   double t69 = -t68 + 2;
   double t70 = -3.0L/2.0L*t17*t4*y1 + (3.0L/2.0L)*y0;
   double t71 = 3*t62 + 3*t64 - t69*t70;
   double t72 = t42*(t50*t60 + t56 + t59) + t48*t67 + t53*t71 + t61*t66;
   double t73 = t1*t7*x1;
   double t74 = t73*z1;
   double t75 = -t74;
   double t76 = t2*t4;
   double t77 = -t76 + 1;
   double t78 = t45*t77;
   double t79 = (1.0L/2.0L)*t4*x1*z1;
   double t80 = -t38*t45 - t39*t63 + t40*t77;
   double t81 = t19*t45;
   double t82 = t21*t63;
   double t83 = -t23*t77 + t81 + t82;
   double t84 = grad*s0*s2*t25*t37*t83;
   double t85 = 2*t76;
   double t86 = -3.0L/2.0L*t17*t4*z1 + (3.0L/2.0L)*z0;
   double t87 = 3*t81 + 3*t82 - t86*(-t85 + 2);
   double t88 = t42*(t50*t79 + t75 + t78) + t48*t80 + t53*t87 + t61*t84;
   double t89 = -1.0L/2.0L*t17*t4*y1 + (1.0L/2.0L)*y0;
   double t90 = t4*x0;
   double t91 = 2*t90;
   double t92 = t17*t7*x1;
   double t93 = 4*t92;
   double t94 = t93*y1;
   double t95 = -t91*y1 + t94;
   double t96 = -1.0L/2.0L*t17*t4*z1 + (1.0L/2.0L)*z0;
   double t97 = t93*z1;
   double t98 = -t91*z1 + t97;
   double t99 = -1.0L/2.0L*t17*t4*x1 + (1.0L/2.0L)*x0;
   double t100 = t14*t4;
   double t101 = t17*t4;
   double t102 = -2*t101;
   double t103 = t0*t17*t7;
   double t104 = -2*t100 + t102 + 4*t103;
   double t105 = -t104*t99 - t89*t95 - t96*t98;
   double t106 = grad*s0*s2*t105*t25*t37;
   double t107 = (1.0L/2.0L)*t29*t4*y1 - 1.0L/2.0L*y2;
   double t108 = t4*x2;
   double t109 = 2*t108;
   double t110 = t29*t7*x1;
   double t111 = 4*t110;
   double t112 = -t111*y1;
   double t113 = t109*y1 + t112;
   double t114 = (1.0L/2.0L)*t29*t4*z1 - 1.0L/2.0L*z2;
   double t115 = -t111*z1;
   double t116 = t109*z1 + t115;
   double t117 = (1.0L/2.0L)*t29*t4*x1 - 1.0L/2.0L*x2;
   double t118 = t26*t4;
   double t119 = t29*t4;
   double t120 = 2*t119;
   double t121 = t0*t29*t7;
   double t122 = 2*t118 + t120 - 4*t121;
   double t123 = -t107*t113 - t114*t116 - t117*t122;
   double t124 = 1.0/(sqrt(t24));
   double t125 = 1.0/(sqrt(t36*t36*t36));
   double t126 = grad*s0*s2*t124*t125*t61;
   double t127 = grad*s0*s2*t124*t37;
   double t128 = -t20*t39;
   double t129 = -t22*t40;
   double t130 = t0*t7*y1;
   double t131 = 2*t130;
   double t132 = t0*t7*z1;
   double t133 = 2*t132;
   double t134 = t108*y1;
   double t135 = 2*t110;
   double t136 = t135*y1;
   double t137 = -t134 + t136;
   double t138 = t108*z1;
   double t139 = t135*z1;
   double t140 = -t138 + t139;
   double t141 = 2*t18;
   double t142 = x1*x1*x1;
   double t143 = t142*t7;
   double t144 = -t141 + 2*t143;
   double t145 = -t119;
   double t146 = 2*t121;
   double t147 = -t118 + t145 + t146;
   double t148 = -t104*t51 - t70*t95 - t86*t98;
   double t149 = grad*s0*s2*t125*t25*t41*t47;
   double t150 = t20*t21;
   double t151 = t22*t23;
   double t152 = 2*t92;
   double t153 = t152*y1;
   double t154 = t153 - t90*y1;
   double t155 = t154*t43;
   double t156 = t152*z1;
   double t157 = t156 - t90*z1;
   double t158 = t157*t45;
   double t159 = 4*t18;
   double t160 = 4*t143;
   double t161 = -t101;
   double t162 = -t100 + 2*t103 + t161;
   double t163 = t137*t21 + t140*t23 + t147*t19 + t154*t39 + t157*t40 + t162*t38;
   double t164 = t106*t61 + t123*t126 + t123*t149 + t127*(t128 + t129 + t131*t39 + t133*t40 - t137*t43 - t140*t45 + t144*t38 + t147*t6) + t148*t53 + t163*t48 + t42*(-t131*t21 - t133*t23 + t150 + t151 + t155 + t158 - t162*t6 - t99*(-t159 + t160));
   double t165 = -t141*y0 + t94;
   double t166 = t4*y0;
   double t167 = 2*t166;
   double t168 = t17*t7*y1;
   double t169 = 4*t168;
   double t170 = t169*z1;
   double t171 = -t167*z1 + t170;
   double t172 = t15*t4;
   double t173 = t1*t17*t7;
   double t174 = t102 - 2*t172 + 4*t173;
   double t175 = -t165*t99 - t171*t96 - t174*t89;
   double t176 = grad*s0*s2*t175*t25*t37;
   double t177 = t112 + t141*y2;
   double t178 = t4*y2;
   double t179 = 2*t178;
   double t180 = t29*t7*y1;
   double t181 = 4*t180;
   double t182 = -t181*z1;
   double t183 = t179*z1 + t182;
   double t184 = t27*t4;
   double t185 = t1*t29*t7;
   double t186 = t120 + 2*t184 - 4*t185;
   double t187 = -t107*t186 - t114*t183 - t117*t177;
   double t188 = 2*t7*x1*y1*z1;
   double t189 = 2*t73;
   double t190 = t131*t38 + t188*t40 + t189*t39;
   double t191 = t178*z1;
   double t192 = 2*t180;
   double t193 = t192*z1;
   double t194 = -t191 + t193;
   double t195 = t18*y2;
   double t196 = t136 - t195;
   double t197 = 2*t185;
   double t198 = t145 - t184 + t197;
   double t199 = -t165*t51 - t171*t86 - t174*t70;
   double t200 = t188*t23;
   double t201 = t189*t21;
   double t202 = t131*t19;
   double t203 = -t200 - t201 - t202;
   double t204 = 2*t168;
   double t205 = t204*z1;
   double t206 = -t166*z1 + t205;
   double t207 = t206*t45;
   double t208 = t153 - t18*y0;
   double t209 = t161 - t172 + 2*t173;
   double t210 = t209*t43;
   double t211 = t19*t196 + t194*t23 + t198*t21 + t206*t40 + t208*t38 + t209*t39;
   double t212 = t126*t187 + t127*(-t18*t39 + t190 - t194*t45 + t196*t6 - t198*t43) + t149*t187 + t176*t61 + t199*t53 + t211*t48 + t42*(t18*t21 + t203 + t207 - t208*t6 + t210);
   double t213 = -t141*z0 + t97;
   double t214 = 2*t20;
   double t215 = t170 - t214*z0;
   double t216 = t16*t4;
   double t217 = t17*t2*t7;
   double t218 = t102 - 2*t216 + 4*t217;
   double t219 = -t213*t99 - t215*t89 - t218*t96;
   double t220 = grad*s0*s2*t219*t25*t37;
   double t221 = t115 + t141*z2;
   double t222 = t182 + t214*z2;
   double t223 = t28*t4;
   double t224 = t2*t29*t7;
   double t225 = t120 + 2*t223 - 4*t224;
   double t226 = -t107*t222 - t114*t225 - t117*t221;
   double t227 = 2*t54;
   double t228 = t133*t38 + t188*t39 + t227*t40;
   double t229 = t20*z2;
   double t230 = t193 - t229;
   double t231 = t18*z2;
   double t232 = t139 - t231;
   double t233 = 2*t224;
   double t234 = t145 - t223 + t233;
   double t235 = -t213*t51 - t215*t70 - t218*t86;
   double t236 = t188*t21;
   double t237 = t227*t23;
   double t238 = t133*t19;
   double t239 = -t236 - t237 - t238;
   double t240 = -t20*z0 + t205;
   double t241 = t240*t43;
   double t242 = t156 - t18*z0;
   double t243 = t161 - t216 + 2*t217;
   double t244 = t243*t45;
   double t245 = t19*t232 + t21*t230 + t23*t234 + t240*t39 + t242*t38 + t243*t40;
   double t246 = t126*t226 + t127*(-t18*t40 + t228 - t230*t43 + t232*t6 - t234*t45) + t149*t226 + t220*t61 + t235*t53 + t245*t48 + t42*(t18*t23 + t239 + t241 - t242*t6 + t244);
   double t247 = t19*t6 - t44 - t46;
   double t248 = grad*s0*s2*t247*t25*t37;
   double t249 = t33*t43;
   double t250 = t35*t45;
   double t251 = t5 - 1;
   double t252 = -t249 - t250 - t251*t31;
   double t253 = t126*t252 + t127*(t11 + t6*t6 + t9) + t149*t252 + t248*t47;
   double t254 = t127*(-t43*t6 + t55 - t59);
   double t255 = t21*t58 - t62 - t64;
   double t256 = grad*s0*s2*t25*t255*t37;
   double t257 = t31*t43;
   double t258 = t35*t63;
   double t259 = t57 - 1;
   double t260 = -t257 - t258 - t259*t33;
   double t261 = grad*s0*s2*t124*t125*t260;
   double t262 = t149*t260 + t254 + t256*t47 + t261*t61;
   double t263 = t127*(-t45*t6 + t74 - t78);
   double t264 = t23*t77 - t81 - t82;
   double t265 = t31*t45;
   double t266 = t33*t63;
   double t267 = t76 - 1;
   double t268 = -t265 - t266 - t267*t35;
   double t269 = t126*t268 + t149*t268 + t263 + t264*t48;
   double t270 = t1*t2*t7;
   double t271 = -t270;
   double t272 = t10 + t271;
   double t273 = grad*s0*s2*t37*t41*t52*t65;
   double t274 = t130*z1;
   double t275 = -t274;
   double t276 = t63*t77;
   double t277 = (1.0L/2.0L)*t4*y1*z1;
   double t278 = t273*t87 + t42*(t275 + t276 + t277*t69) + t66*t80 + t67*t84;
   double t279 = grad*s0*s2*t124*t125*t67;
   double t280 = grad*s0*s2*t125*t25*t41*t65;
   double t281 = t157*t63;
   double t282 = t162*t43;
   double t283 = t106*t67 + t123*t279 + t123*t280 + t127*(t137*t58 - t140*t63 - t147*t43 + t190 - t20*t38) + t148*t273 + t163*t66 + t42*(-t154*t58 + t19*t20 + t203 + t281 + t282);
   double t284 = -t18*t38;
   double t285 = t1*t7*z1;
   double t286 = 2*t285;
   double t287 = y1*y1*y1;
   double t288 = t287*t7;
   double t289 = -t214 + 2*t288;
   double t290 = t18*t19;
   double t291 = t208*t43;
   double t292 = t206*t63;
   double t293 = 4*t20;
   double t294 = 4*t288;
   double t295 = t127*(t129 + t189*t38 - t194*t63 - t196*t43 + t198*t58 + t284 + t286*t40 + t289*t39) + t176*t67 + t187*t279 + t187*t280 + t199*t273 + t211*t66 + t42*(t151 - t189*t19 - t209*t58 - t23*t286 + t290 + t291 + t292 - t89*(-t293 + t294));
   double t296 = t2*t7*y1;
   double t297 = 2*t296;
   double t298 = t188*t38 + t286*t39 + t297*t40;
   double t299 = t188*t19;
   double t300 = t23*t297;
   double t301 = t21*t286;
   double t302 = -t299 - t300 - t301;
   double t303 = t242*t43;
   double t304 = t243*t63;
   double t305 = t127*(-t20*t40 + t230*t58 - t232*t43 - t234*t63 + t298) + t220*t67 + t226*t279 + t226*t280 + t235*t273 + t245*t66 + t42*(t20*t23 - t240*t58 + t302 + t303 + t304);
   double t306 = t247*t66 + t252*t279 + t252*t280 + t254;
   double t307 = t127*(t270 + t58*t58 + t9) + t256*t65 + t260*t279 + t260*t280;
   double t308 = t127*(t274 - t276 - t58*t63);
   double t309 = t264*t66 + t268*t279 + t268*t280 + t308;
   double t310 = t12 + t271;
   double t311 = grad*s0*s2*t37*t41*t52*t83;
   double t312 = grad*s0*s2*t123*t124*t125;
   double t313 = grad*s0*s2*t125*t25*t41*t83;
   double t314 = t154*t63;
   double t315 = t162*t45;
   double t316 = t106*t80 + t123*t313 + t127*(-t137*t63 + t140*t77 - t147*t45 - t22*t38 + t228) + t148*t311 + t163*t84 + t312*t80 + t42*(-t157*t77 + t19*t22 + t239 + t314 + t315);
   double t317 = grad*s0*s2*t124*t125*t187;
   double t318 = t208*t45;
   double t319 = t209*t63;
   double t320 = t127*(t194*t77 - t196*t45 - t198*t63 - t22*t39 + t298) + t176*t80 + t187*t313 + t199*t311 + t211*t84 + t317*t80 + t42*(-t206*t77 + t21*t22 + t302 + t318 + t319);
   double t321 = grad*s0*s2*t124*t125*t226;
   double t322 = 2*t22;
   double t323 = z1*z1*z1;
   double t324 = t323*t7;
   double t325 = -t322 + 2*t324;
   double t326 = t242*t45;
   double t327 = t240*t63;
   double t328 = 4*t22;
   double t329 = 4*t324;
   double t330 = t127*(t128 + t227*t38 - t230*t63 - t232*t45 + t234*t77 + t284 + t297*t39 + t325*t40) + t220*t80 + t226*t313 + t235*t311 + t245*t84 + t321*t80 + t42*(t150 - t19*t227 - t21*t297 - t243*t77 + t290 + t326 + t327 - t96*(-t328 + t329));
   double t331 = grad*s0*s2*t124*t125*t252;
   double t332 = t248*t83 + t252*t313 + t263 + t331*t80;
   double t333 = t256*t83 + t260*t313 + t261*t80 + t308;
   double t334 = grad*s0*s2*t124*t125*t268;
   double t335 = t127*(t11 + t270 + t77*t77) + t264*t84 + t268*t313 + t334*t80;
   double t336 = grad*s0*s2*t105*t125*t25*t41;
   double t337 = 1.0/(sqrt(t36*t36*t36*t36*t36));
   double t338 = (3.0L/2.0L)*t29*t4*y1 - 3.0L/2.0L*y2;
   double t339 = (3.0L/2.0L)*t29*t4*z1 - 3.0L/2.0L*z2;
   double t340 = (3.0L/2.0L)*t29*t4*x1 - 3.0L/2.0L*x2;
   double t341 = 2*t137*t21 + 2*t140*t23 + 2*t147*t19 + 2*t154*t39 + 2*t157*t40 + 2*t162*t38;
   double t342 = grad*s0*s2*t124*t125*t41;
   double t343 = -t136;
   double t344 = t134 + t343;
   double t345 = -t139;
   double t346 = t138 + t345;
   double t347 = t7*x1*x2*y1;
   double t348 = -t181;
   double t349 = 1.0/(t3*t3*t3);
   double t350 = t0*t29*t349*y1;
   double t351 = t348 + 16*t350;
   double t352 = t7*x1*x2*z1;
   double t353 = t29*t7*z1;
   double t354 = -4*t353;
   double t355 = t0*t29*t349*z1;
   double t356 = t354 + 16*t355;
   double t357 = t0*t7*x2;
   double t358 = t142*t29*t349;
   double t359 = t118 + t119 - t146;
   double t360 = t7*x0*x1*y1;
   double t361 = t0*t17*t349*y1;
   double t362 = t169 - 16*t361;
   double t363 = t7*x0*x1*z1;
   double t364 = t17*t7*z1;
   double t365 = 4*t364;
   double t366 = t0*t17*t349*z1;
   double t367 = t365 - 16*t366;
   double t368 = t0*t7*x0;
   double t369 = t142*t17*t349;
   double t370 = 4*t347;
   double t371 = t192 - 8*t350;
   double t372 = 4*t352;
   double t373 = 2*t353;
   double t374 = -8*t355 + t373;
   double t375 = 4*t360;
   double t376 = t204 - 8*t361;
   double t377 = 4*t363;
   double t378 = 2*t364;
   double t379 = -8*t366 + t378;
   double t380 = -t91;
   double t381 = grad*s0*s2*t199*t37*t41*t52;
   double t382 = grad*s0*s2*t125*t175*t25*t41;
   double t383 = grad*s0*s2*t124*t337*t41*(-t177*t340 - t183*t339 - t186*t338);
   double t384 = -t193;
   double t385 = t191 + t384;
   double t386 = t7*x1*y2*z1;
   double t387 = -4*t386;
   double t388 = t7*x2*y1*z1;
   double t389 = -4*t388;
   double t390 = t29*t349*x1*y1*z1;
   double t391 = 16*t390;
   double t392 = t119 + t184 - t197;
   double t393 = t195 + t343;
   double t394 = t0*t7*y2;
   double t395 = t7*x1*y1*y2;
   double t396 = 4*t395;
   double t397 = t1*t7*x2;
   double t398 = -t111;
   double t399 = t1*t29*t349*x1;
   double t400 = t398 + 16*t399;
   double t401 = t7*x1*y0*z1;
   double t402 = 4*t401;
   double t403 = t7*x0*y1*z1;
   double t404 = t17*t349*x1*y1*z1;
   double t405 = -16*t404;
   double t406 = 4*t403 + t405;
   double t407 = -t167;
   double t408 = t0*t7*y0;
   double t409 = t7*x1*y0*y1;
   double t410 = 4*t409;
   double t411 = t1*t7*x0;
   double t412 = t1*t17*t349*x1;
   double t413 = -16*t412 + t93;
   double t414 = 2*t386;
   double t415 = 2*t388;
   double t416 = -8*t390;
   double t417 = 2*t401;
   double t418 = -8*t404;
   double t419 = 2*t403 + t418;
   double t420 = -t178;
   double t421 = -t108;
   double t422 = t135 - 8*t399;
   double t423 = -t166;
   double t424 = -t90;
   double t425 = t152 - 8*t412;
   double t426 = t105*t381 + t106*t211 + t123*t382 + t123*t383 + t127*(t137*t209 + t140*t206 + t147*t208 + t154*t198 + t157*t194 + t162*t196 + t19*(2*t347 + t371 + 2*t394 + t420) + t21*(2*t395 + 2*t397 + t421 + t422) + t23*(t414 + t415 + t416) + t38*(2*t360 + t376 + 2*t408 + t423) + t39*(2*t409 + 2*t411 + t424 + t425) + t40*(t417 + t419)) + t163*t176 + t163*t317 + t187*t336 + t211*t312 + t342*(-t107*(t109 - t396 - 4*t397 + t400) - t114*(t387 + t389 + t391) - t117*(t179 + t351 - t370 - 4*t394) - t344*t392 - t346*t385 - t359*t393) + t42*(-t154*t209 - t157*t206 - t162*t208 - t89*(t380 + t410 + 4*t411 + t413) - t96*(t402 + t406) - t99*(t362 + t375 + t407 + 4*t408));
   double t427 = grad*s0*s2*t235*t37*t41*t52;
   double t428 = grad*s0*s2*t125*t219*t25*t41;
   double t429 = grad*s0*s2*t124*t337*t41*(-t221*t340 - t222*t338 - t225*t339);
   double t430 = t7*x1*y1*z2;
   double t431 = t391 - 4*t430;
   double t432 = t119 + t223 - t233;
   double t433 = -t29*t7*x1*z1 + (1.0L/2.0L)*t4*x1*z2;
   double t434 = t4*z2;
   double t435 = 2*t434;
   double t436 = t0*t7*z2;
   double t437 = t7*x1*z1*z2;
   double t438 = 4*t437;
   double t439 = t2*t7*x2;
   double t440 = t2*t29*t349*x1;
   double t441 = t398 + 16*t440;
   double t442 = t7*x1*y1*z0;
   double t443 = 4*t442;
   double t444 = t17*t7*x1*z1 - 1.0L/2.0L*t4*x1*z0;
   double t445 = t4*z0;
   double t446 = -2*t445;
   double t447 = t0*t7*z0;
   double t448 = t7*x1*z0*z1;
   double t449 = 4*t448;
   double t450 = t2*t7*x0;
   double t451 = t17*t2*t349*x1;
   double t452 = -16*t451 + t93;
   double t453 = t416 + 2*t430;
   double t454 = 2*t442;
   double t455 = -t434;
   double t456 = t135 - 8*t440;
   double t457 = -t445;
   double t458 = t152 - 8*t451;
   double t459 = t105*t427 + t106*t245 + t123*t428 + t123*t429 + t127*(t137*t240 + t140*t243 + t147*t242 + t154*t230 + t157*t234 + t162*t232 + t19*(2*t352 + t374 + 2*t436 + t455) + t21*(t415 + t453) + t23*(t421 + 2*t437 + 2*t439 + t456) + t38*(2*t363 + t379 + 2*t447 + t457) + t39*(t419 + t454) + t40*(t424 + 2*t448 + 2*t450 + t458)) + t163*t220 + t163*t321 + t226*t336 + t245*t312 + t342*(-t107*(t389 + t431) - t114*(t109 - t438 - 4*t439 + t441) - t117*(t356 - t372 + t435 - 4*t436) - t122*t433 - t344*(t229 + t384) - t346*t432) + t42*(-t104*t444 - t154*t240 - t157*t243 - t89*(t406 + t443) - t96*(t380 + t449 + 4*t450 + t452) - t99*(t367 + t377 + t446 + 4*t447));
   double t460 = grad*s0*s2*t124*t125*t247;
   double t461 = -t20;
   double t462 = t131 + t461;
   double t463 = -t22;
   double t464 = t133 + t463;
   double t465 = t49 - 2;
   double t466 = grad*s0*s2*t124*t337*t41*(-3*t249 - 3*t250 - t340*t465);
   double t467 = -4*t130 + t214;
   double t468 = -4*t132 + t322;
   double t469 = (1.0L/2.0L)*t0*t4 - 1.0L/2.0L;
   double t470 = t105*t248 + t123*t460 + t123*t466 + t127*(t144*t19 - t155 - t158 + t162*t6 + t21*t462 + t23*t464) + t163*t331 + t252*t336 + t342*(-t107*t467 - t113*t60 - t114*t468 - t116*t79 - t117*(t159 - t160) - t122*t469);
   double t471 = grad*s0*s2*t124*t125*t255;
   double t472 = t188*t35;
   double t473 = (1.0L/2.0L)*t1*t4 - 1.0L/2.0L;
   double t474 = t68 - 2;
   double t475 = grad*s0*s2*t124*t337*t41*(-3*t257 - 3*t258 - t338*t474);
   double t476 = t105*t256 + t123*t471 + t123*t475 + t127*(t154*t58 + t19*t462 + t200 + t201 - t281 - t282) + t163*t261 + t260*t336 + t342*(-t113*t473 - t116*t277 - t117*t467 - t122*t60 + t189*t33 + t472);
   double t477 = t188*t33;
   double t478 = (1.0L/2.0L)*t2*t4 - 1.0L/2.0L;
   double t479 = grad*s0*s2*t124*t337*t41*(-3*t265 - 3*t266 - t339*(t85 - 2));
   double t480 = t106*t264 + t123*t479 + t127*(t157*t77 + t19*t464 + t236 + t237 - t314 - t315) + t163*t334 + t264*t312 + t268*t336 + t342*(-t113*t277 - t116*t478 - t117*t468 - t122*t79 + t227*t35 + t477);
   double t481 = 2*t19*t196 + 2*t194*t23 + 2*t198*t21 + 2*t206*t40 + 2*t208*t38 + 2*t209*t39;
   double t482 = t7*y1*y2*z1;
   double t483 = t1*t29*t349*z1;
   double t484 = t354 + 16*t483;
   double t485 = t1*t7*y2;
   double t486 = t287*t29*t349;
   double t487 = t7*y0*y1*z1;
   double t488 = t1*t17*t349*z1;
   double t489 = t365 - 16*t488;
   double t490 = t1*t7*y0;
   double t491 = t17*t287*t349;
   double t492 = 4*t482;
   double t493 = t373 - 8*t483;
   double t494 = 4*t487;
   double t495 = t378 - 8*t488;
   double t496 = -t29*t7*y1*z1 + (1.0L/2.0L)*t4*y1*z2;
   double t497 = t1*t7*z2;
   double t498 = t7*y1*z1*z2;
   double t499 = 4*t498;
   double t500 = t2*t7*y2;
   double t501 = t2*t29*t349*y1;
   double t502 = t348 + 16*t501;
   double t503 = t17*t7*y1*z1 - 1.0L/2.0L*t4*y1*z0;
   double t504 = t1*t7*z0;
   double t505 = t7*y1*z0*z1;
   double t506 = 4*t505;
   double t507 = t2*t7*y0;
   double t508 = t17*t2*t349*y1;
   double t509 = t169 - 16*t508;
   double t510 = t192 - 8*t501;
   double t511 = t204 - 8*t508;
   double t512 = t127*(t19*(t414 + t453) + t194*t243 + t196*t242 + t198*t240 + t206*t234 + t208*t232 + t209*t230 + t21*(t455 + 2*t482 + t493 + 2*t497) + t23*(t420 + 2*t498 + 2*t500 + t510) + t38*(t417 + t418 + t454) + t39*(t457 + 2*t487 + t495 + 2*t504) + t40*(t423 + 2*t505 + 2*t507 + t511)) + t175*t427 + t176*t245 + t187*t428 + t187*t429 + t211*t220 + t211*t321 + t226*t382 + t245*t317 + t342*(-t107*(t435 + t484 - t492 - 4*t497) - t114*(t179 - t499 - 4*t500 + t502) - t117*(t387 + t431) - t186*t496 - t385*t432 - t393*(t231 + t345)) + t42*(-t174*t503 - t206*t243 - t208*t242 - t89*(t446 + t489 + t494 + 4*t504) - t96*(t407 + t506 + 4*t507 + t509) - t99*(t402 + t405 + t443));
   double t513 = -t18;
   double t514 = t189 + t513;
   double t515 = t141 - 4*t73;
   double t516 = t127*(t200 + t202 - t207 + t208*t6 + t21*t514 - t210) + t175*t248 + t187*t460 + t187*t466 + t211*t331 + t252*t382 + t342*(-t107*t515 + t131*t31 - t177*t469 - t183*t79 - t186*t60 + t472);
   double t517 = t286 + t463;
   double t518 = -4*t285 + t322;
   double t519 = t127*(t19*t514 + t209*t58 + t21*t289 + t23*t517 - t291 - t292) + t175*t256 + t187*t471 + t187*t475 + t211*t261 + t260*t382 + t342*(-t107*(t293 - t294) - t114*t518 - t117*t515 - t177*t60 - t183*t277 - t186*t473);
   double t520 = t188*t31;
   double t521 = t127*(t206*t77 + t21*t517 + t299 + t300 - t318 - t319) + t176*t264 + t187*t479 + t211*t334 + t264*t317 + t268*t382 + t342*(-t107*t518 - t177*t79 - t183*t478 - t186*t277 + t297*t35 + t520);
   double t522 = 2*t19*t232 + 2*t21*t230 + 2*t23*t234 + 2*t240*t39 + 2*t242*t38 + 2*t243*t40;
   double t523 = t2*t7*z2;
   double t524 = t29*t323*t349;
   double t525 = t2*t7*z0;
   double t526 = t17*t323*t349;
   double t527 = t227 + t513;
   double t528 = t141 - 4*t54;
   double t529 = t127*(t23*t527 + t236 + t238 - t241 + t242*t6 - t244) + t219*t248 + t226*t460 + t226*t466 + t245*t331 + t252*t428 + t342*(-t114*t528 + t133*t31 - t221*t469 - t222*t60 - t225*t79 + t477);
   double t530 = t297 + t461;
   double t531 = t214 - 4*t296;
   double t532 = t127*(t23*t530 + t240*t58 + t299 + t301 - t303 - t304) + t219*t256 + t226*t471 + t226*t475 + t245*t261 + t260*t428 + t342*(-t114*t531 - t221*t60 - t222*t473 - t225*t277 + t286*t33 + t520);
   double t533 = t127*(t19*t527 + t21*t530 + t23*t325 + t243*t77 - t326 - t327) + t220*t264 + t226*t479 + t245*t334 + t264*t321 + t268*t428 + t342*(-t107*t531 - t114*(t328 - t329) - t117*t528 - t221*t79 - t222*t277 - t225*t478);
   double t534 = t252*t471 + t252*t475 + t260*t460 + t342*(-t259*t43 - t465*t60 + t56);
   double t535 = t267*z1;
   double t536 = t252*t479 + t264*t331 + t268*t460 + t342*(-t18*t535 - t465*t79 + t75);
   double t537 = t260*t479 + t261*t264 + t268*t471 + t342*(-t20*t535 + t275 - t277*t474);
   // Fill up Hessian
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i0+0] +=  t42*(t13 - pow(t6, 2)) + t48*(2*t38*t6 - 2*t39*t4*x1*y1 - 2*t4*t40*x1*z1) + t53*(3*t44 + 3*t46 - t50*t51);
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i0+1] +=  t72;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i0+2] +=  t88;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i0+0] +=  t72;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i0+1] +=  t273*t71 + t42*(t272 - pow(t58, 2)) + t66*(-2*t38*t4*x1*y1 + 2*t39*t58 - 2*t4*t40*y1*z1);
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i0+2] +=  t278;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i0+0] +=  t88;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i0+1] +=  t278;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i0+2] +=  t311*t87 + t42*(t310 - pow(t77, 2)) + t84*(-2*t38*t4*x1*z1 - 2*t39*t4*y1*z1 + 2*t40*t77);
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i1+0] +=  t164;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i1+1] +=  t212;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i1+2] +=  t246;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i1+0] +=  t283;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i1+1] +=  t295;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i1+2] +=  t305;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i1+0] +=  t316;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i1+1] +=  t320;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i1+2] +=  t330;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i2+0] +=  t253;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i2+1] +=  t262;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i2+2] +=  t269;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i2+0] +=  t306;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i2+1] +=  t307;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i2+2] +=  t309;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i2+0] +=  t332;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i2+1] +=  t333;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i2+2] +=  t335;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i0+0] +=  t164;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i0+1] +=  t283;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i0+2] +=  t316;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i0+0] +=  t212;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i0+1] +=  t295;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i0+2] +=  t320;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i0+0] +=  t246;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i0+1] +=  t305;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i0+2] +=  t330;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i1+0] +=  grad*s0*s2*t105*t148*t37*t41*t52 + grad*s0*s2*t123*t124*t337*t41*(-t113*t338 - t116*t339 - t122*t340) + t106*t341 + t127*(2*t137*t154 + 2*t140*t157 + 2*t147*t162 + t19*(-t109 + 6*t110 + 4*t357 - 8*t358) + t21*(t370 + t371) + t23*(t372 + t374) + t38*(4*t368 - 8*t369 + t380 + 6*t92) + t39*(t375 + t376) + t40*(t377 + t379)) + t312*t341 + t336*(-t113*t33 - t116*t35 - t122*t31) + t342*(-t107*(-8*t347 + t351) - t114*(-8*t352 + t356) - t117*(4*t108 - 12*t110 - 8*t357 + 16*t358) - pow(t344, 2) - pow(t346, 2) - pow(t359, 2)) + t42*(-pow(t154, 2) - pow(t157, 2) - pow(t162, 2) - t89*(8*t360 + t362) - t96*(8*t363 + t367) - t99*(8*t368 - 16*t369 - 4*t90 + 12*t92));
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i1+1] +=  t426;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i1+2] +=  t459;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i1+0] +=  t426;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i1+1] +=  t127*(t19*(t396 + t422) + 2*t194*t206 + 2*t196*t208 + 2*t198*t209 + t21*(-t179 + 6*t180 + 4*t485 - 8*t486) + t23*(t492 + t493) + t38*(t410 + t425) + t39*(6*t168 + t407 + 4*t490 - 8*t491) + t40*(t494 + t495)) + t175*t381 + t176*t481 + t187*t383 + t317*t481 + t342*(-t107*(4*t178 - 12*t180 - 8*t485 + 16*t486) - t114*(-8*t482 + t484) - t117*(-8*t395 + t400) - 1.0L/2.0L*t186*t392 - pow(t385, 2) - pow(t393, 2)) + t382*(-t177*t31 - t183*t35 - t186*t33) + t42*(-1.0L/2.0L*t174*t209 - pow(t206, 2) - pow(t208, 2) - t89*(-4*t166 + 12*t168 + 8*t490 - 16*t491) - t96*(8*t487 + t489) - t99*(8*t409 + t413));
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i1+2] +=  t512;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i1+0] +=  t459;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i1+1] +=  t512;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i1+2] +=  t127*(t19*(t438 + t456) + t21*(t499 + t510) + t23*(6*t353 - t435 + 4*t523 - 8*t524) + 2*t230*t240 + 2*t232*t242 + 2*t234*t243 + t38*(t449 + t458) + t39*(t506 + t511) + t40*(6*t364 + t446 + 4*t525 - 8*t526)) + t219*t427 + t220*t522 + t226*t429 + t321*t522 + t342*(-t107*(-8*t498 + t502) - t114*(-12*t353 + 4*t434 - 8*t523 + 16*t524) - t117*(-8*t437 + t441) - t221*t433 - t222*t496 - 1.0L/2.0L*t225*t432) + t42*(-t213*t444 - t215*t503 - 1.0L/2.0L*t218*t243 - t89*(8*t505 + t509) - t96*(12*t364 - 4*t445 + 8*t525 - 16*t526) - t99*(8*t448 + t452)) + t428*(-t221*t31 - t222*t33 - t225*t35);
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i2+0] +=  t470;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i2+1] +=  t476;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i2+2] +=  t480;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i2+0] +=  t516;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i2+1] +=  t519;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i2+2] +=  t521;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i2+0] +=  t529;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i2+1] +=  t532;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i2+2] +=  t533;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i0+0] +=  t253;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i0+1] +=  t306;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i0+2] +=  t332;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i0+0] +=  t262;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i0+1] +=  t307;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i0+2] +=  t333;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i0+0] +=  t269;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i0+1] +=  t309;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i0+2] +=  t335;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i1+0] +=  t470;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i1+1] +=  t516;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i1+2] +=  t529;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i1+0] +=  t476;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i1+1] +=  t519;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i1+2] +=  t532;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i1+0] +=  t480;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i1+1] +=  t521;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i1+2] +=  t533;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i2+0] +=  t252*t466 + t342*(t13 - pow(t251, 2)) + t460*(-t31*t465 - 2*t33*t4*x1*y1 - 2*t35*t4*x1*z1);
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i2+1] +=  t534;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i2+2] +=  t536;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i2+0] +=  t534;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i2+1] +=  t260*t475 + t342*(-pow(t259, 2) + t272) + t471*(-2*t31*t4*x1*y1 - t33*t474 - 2*t35*t4*y1*z1);
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i2+2] +=  t537;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i2+0] +=  t536;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i2+1] +=  t537;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i2+2] +=  t268*t479 + t334*(-2*t19*t4*x1*z1 - 2*t21*t4*y1*z1 + 2*t23*t77) + t342*(-pow(t267, 2) + t310);
}

void hessian_dihed_angle(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double* hessian, long ndelta) {
   dlist_row_type *delta0, *delta1, *delta2;
   double tmp = sin(value);
   if (tmp!=0.0) tmp = -grad/tmp;
   // First part
   hessian_dihed_cos(ic, deltas, cos(value), tmp, hessian, ndelta);
   // Second part
   double jacobian[9];
   delta0 = deltas + (*ic).i0;
   delta1 = deltas + (*ic).i1;
   delta2 = deltas + (*ic).i2;
   (*delta0).gx = 0.0;
   (*delta0).gy = 0.0;
   (*delta0).gz = 0.0;
   (*delta1).gx = 0.0;
   (*delta1).gy = 0.0;
   (*delta1).gz = 0.0;
   (*delta2).gx = 0.0;
   (*delta2).gy = 0.0;
   (*delta2).gz = 0.0;
   ic_back_fns[(*ic).kind](ic, deltas, value, 1.0);
   jacobian[0] = (*delta0).gx;
   jacobian[1] = (*delta0).gy;
   jacobian[2] = (*delta0).gz;
   jacobian[3] = (*delta1).gx;
   jacobian[4] = (*delta1).gy;
   jacobian[5] = (*delta1).gz;
   jacobian[6] = (*delta2).gx;
   jacobian[7] = (*delta2).gy;
   jacobian[8] = (*delta2).gz;
   long i,j;
   tmp *= -cos(value);
   for (i=0;i<3;i++){
     for (j=0;j<3;j++){
       hessian[(3*(*ic).i0+i)*ndelta + 3*(*ic).i0+j] -= tmp*jacobian[i]*jacobian[j];
       hessian[(3*(*ic).i0+i)*ndelta + 3*(*ic).i1+j] -= tmp*jacobian[i]*jacobian[j+3];
       hessian[(3*(*ic).i0+i)*ndelta + 3*(*ic).i2+j] -= tmp*jacobian[i]*jacobian[j+6];
       hessian[(3*(*ic).i1+i)*ndelta + 3*(*ic).i0+j] -= tmp*jacobian[i+3]*jacobian[j];
       hessian[(3*(*ic).i1+i)*ndelta + 3*(*ic).i1+j] -= tmp*jacobian[i+3]*jacobian[j+3];
       hessian[(3*(*ic).i1+i)*ndelta + 3*(*ic).i2+j] -= tmp*jacobian[i+3]*jacobian[j+6];
       hessian[(3*(*ic).i2+i)*ndelta + 3*(*ic).i0+j] -= tmp*jacobian[i+6]*jacobian[j];
       hessian[(3*(*ic).i2+i)*ndelta + 3*(*ic).i1+j] -= tmp*jacobian[i+6]*jacobian[j+3];
       hessian[(3*(*ic).i2+i)*ndelta + 3*(*ic).i2+j] -= tmp*jacobian[i+6]*jacobian[j+6];
     }
   }
}

void hessian_oop_cos_low(iclist_row_type* ic, long i0, long i1, long i2, dlist_row_type* delta0, dlist_row_type* delta1, dlist_row_type* delta2, double value, double grad, double* hessian, long ndelta) {
   // Initialization
   double x0 = (*delta0).dx;
   double y0 = (*delta0).dy;
   double z0 = (*delta0).dz;
   double x1 = (*delta1).dx;
   double y1 = (*delta1).dy;
   double z1 = (*delta1).dz;
   double x2 = (*delta2).dx;
   double y2 = (*delta2).dy;
   double z2 = (*delta2).dz;
   // Compute auxiliary variables
   double t0 = y0*z1;
   double t1 = y1*z0;
   double t2 = t0 - t1;
   double t3 = x1*z0;
   double t4 = x0*z1;
   double t5 = t3 - t4;
   double t6 = x0*y1;
   double t7 = x1*y0;
   double t8 = t6 - t7;
   double t9 = t2*x2 + t5*y2 + t8*z2;
   double t10 = y1*z2;
   double t11 = y2*z1;
   double t12 = 2*t10 - 2*t11;
   double t13 = x2*x2;
   double t14 = y2*y2;
   double t15 = z2*z2;
   double t16 = t13 + t14 + t15;
   double t17 = 1.0/t16;
   double t18 = -t3 + t4;
   double t19 = t18*t18 + t2*t2 + t8*t8;
   double t20 = 1.0/t19;
   double t21 = (1.0L/2.0L)*t12*t17*t20;
   double t22 = t21*t9;
   double t23 = t8*y1;
   double t24 = t18*z1;
   double t25 = -2*t23 - 2*t24;
   double t26 = t9*t9;
   double t27 = 1.0/(t19*t19);
   double t28 = (1.0L/2.0L)*t17*t26*t27;
   double t29 = t25*t28;
   double t30 = t20*t26;
   double t31 = -t17*t30 + 1;
   double t32 = 1.0/(sqrt(t31*t31*t31));
   double t33 = -t22 - t29;
   double t34 = grad*t32*t33;
   double t35 = grad/sqrt(t31);
   double t36 = t10 - t11;
   double t37 = -2*pow(y1, 2);
   double t38 = -2*pow(z1, 2);
   double t39 = t17*t25*t27*t9;
   double t40 = 1.0/(t19*t19*t19);
   double t41 = (1.0L/2.0L)*t17*t25*t26*t40;
   double t42 = x1*z2;
   double t43 = x2*z1;
   double t44 = -2*t42 + 2*t43;
   double t45 = (1.0L/2.0L)*t17*t20*t9;
   double t46 = t44*t45;
   double t47 = t8*x1;
   double t48 = t2*z1;
   double t49 = 2*t47 - 2*t48;
   double t50 = t28*t49;
   double t51 = t46 + t50;
   double t52 = -t42 + t43;
   double t53 = (1.0L/2.0L)*t17*t20*t52;
   double t54 = t17*t26*t27*x1;
   double t55 = (1.0L/2.0L)*t17*t27*t49*t9;
   double t56 = 4*t47 - 4*t48;
   double t57 = t34*t51 + t35*(-t12*t53 - t12*t55 - t39*t52 - t41*t56 - t54*y1);
   double t58 = x1*y2;
   double t59 = x2*y1;
   double t60 = 2*t58 - 2*t59;
   double t61 = t45*t60;
   double t62 = t18*x1;
   double t63 = t2*y1;
   double t64 = 2*t62 + 2*t63;
   double t65 = t28*t64;
   double t66 = t61 + t65;
   double t67 = grad*t32*t66;
   double t68 = t58 - t59;
   double t69 = (1.0L/2.0L)*t17*t20*t68;
   double t70 = (1.0L/2.0L)*t17*t27*t64*t9;
   double t71 = 4*t62 + 4*t63;
   double t72 = t33*t67 + t35*(-t12*t69 - t12*t70 - t39*t68 - t41*t71 - t54*z1);
   double t73 = y0*z2;
   double t74 = y2*z0;
   double t75 = -2*t73 + 2*t74;
   double t76 = t45*t75;
   double t77 = 2*y0;
   double t78 = 2*z0;
   double t79 = t18*t78 + t77*t8;
   double t80 = t28*t79;
   double t81 = t76 + t80;
   double t82 = -t73 + t74;
   double t83 = (1.0L/2.0L)*t17*t20*t82;
   double t84 = t77*y1;
   double t85 = t78*z1;
   double t86 = t17*t27*t79*t9;
   double t87 = 4*y0;
   double t88 = 4*z0;
   double t89 = t18*t88 + t8*t87;
   double t90 = t34*t81 + t35*(-t12*t83 - t28*(t84 + t85) - t36*t86 - t39*t82 - t41*t89);
   double t91 = x0*z2;
   double t92 = x2*z0;
   double t93 = 2*t91 - 2*t92;
   double t94 = t45*t93;
   double t95 = 2*x0;
   double t96 = t2*t78 - t8*t95;
   double t97 = t28*t96;
   double t98 = grad*t32*(t94 + t97);
   double t99 = t91 - t92;
   double t100 = (1.0L/2.0L)*t17*t20*t99;
   double t101 = t17*t20*t9;
   double t102 = t101*z2;
   double t103 = 4*t6;
   double t104 = 2*t7;
   double t105 = (1.0L/2.0L)*t17*t27*t9*t96;
   double t106 = 4*x0;
   double t107 = (1.0L/2.0L)*t17*t26*t40*(-t106*t8 + t2*t88);
   double t108 = t33*t98 + t35*(-t100*t12 - t102 - t105*t12 - t107*t25 - t28*(-t103 + t104) - t39*t99);
   double t109 = x0*y2;
   double t110 = x2*y0;
   double t111 = -2*t109 + 2*t110;
   double t112 = (1.0L/2.0L)*t111*t17*t20;
   double t113 = t112*t9;
   double t114 = -t18*t95 - t2*t77;
   double t115 = t114*t28;
   double t116 = grad*t32*(t113 + t115);
   double t117 = -t109 + t110;
   double t118 = (1.0L/2.0L)*t117*t17*t20;
   double t119 = t101*y2;
   double t120 = -4*t4;
   double t121 = 2*t3;
   double t122 = (1.0L/2.0L)*t114*t17*t27*t9;
   double t123 = (1.0L/2.0L)*t17*t26*t40*(-t106*t18 - t2*t87);
   double t124 = t116*t33 + t35*(-t117*t39 - t118*t12 + t119 - t12*t122 - t123*t25 - t28*(t120 + t121));
   double t125 = 1.0/(t16*t16);
   double t126 = t125*t30;
   double t127 = t126*x2;
   double t128 = 2*t0;
   double t129 = 2*t1;
   double t130 = t128 - t129;
   double t131 = t130*t45;
   double t132 = -t127 + t131;
   double t133 = grad*t132*t32;
   double t134 = (1.0L/2.0L)*t17*t2*t20;
   double t135 = t125*t20*t9*x2;
   double t136 = t125*t26*t27*x2;
   double t137 = t133*t33 + t35*(-t12*t134 + t12*t135 + t136*t25 - t2*t39);
   double t138 = t126*y2;
   double t139 = 2*t4;
   double t140 = t121 - t139;
   double t141 = t140*t45;
   double t142 = -t138 + t141;
   double t143 = grad*t142*t32;
   double t144 = (1.0L/2.0L)*t17*t20*t5;
   double t145 = t101*z1;
   double t146 = t125*t20*t9*y2;
   double t147 = t125*t26*t27*y2;
   double t148 = t143*t33 + t35*(-t12*t144 + t12*t146 + t145 + t147*t25 - t39*t5);
   double t149 = t126*z2;
   double t150 = 2*t6;
   double t151 = -t104 + t150;
   double t152 = t151*t45;
   double t153 = -t149 + t152;
   double t154 = grad*t153*t32;
   double t155 = (1.0L/2.0L)*t17*t20*t8;
   double t156 = t101*y1;
   double t157 = t125*t20*t9*z2;
   double t158 = t125*t26*t27*z2;
   double t159 = t154*t33 + t35*(-t12*t155 + t12*t157 - t156 + t158*t25 - t39*t8);
   double t160 = -t46 - t50;
   double t161 = grad*t160*t32;
   double t162 = -2*pow(x1, 2);
   double t163 = t17*t27*t49*t9;
   double t164 = (1.0L/2.0L)*t17*t26*t40*t49;
   double t165 = t17*t26*t27;
   double t166 = t17*t27*t64*t9;
   double t167 = t160*t67 + t35*(-t163*t68 - t164*t71 - t165*y1*z1 - t166*t52 - t44*t69);
   double t168 = -4*t7;
   double t169 = t161*t81 + t35*(t102 - t164*t89 - t28*(t150 + t168) - t44*t83 - t52*t86 - t55*t75);
   double t170 = t95*x1;
   double t171 = t160*t98 + t35*(-t100*t44 - t105*t44 - t107*t49 - t163*t99 - t28*(t170 + t85));
   double t172 = t101*x2;
   double t173 = 4*t0;
   double t174 = t116*t160 + t35*(-t117*t163 - t118*t44 - t122*t44 - t123*t49 - t172 - t28*(t129 - t173));
   double t175 = t133*t160 + t35*(-t130*t55 - t134*t44 + t135*t44 + t136*t49 - t145);
   double t176 = t143*t160 + t35*(-t144*t44 + t146*t44 + t147*t49 - t163*t5);
   double t177 = t101*x1;
   double t178 = t154*t160 + t35*(-t155*t44 + t157*t44 + t158*t49 - t163*t8 + t177);
   double t179 = -t61 - t65;
   double t180 = grad*t179*t32;
   double t181 = (1.0L/2.0L)*t17*t26*t40*t64;
   double t182 = 4*t3;
   double t183 = t180*t81 + t35*(-t119 - t181*t89 - t28*(t139 - t182) - t60*t83 - t68*t86 - t70*t75);
   double t184 = -4*t1;
   double t185 = t179*t98 + t35*(-t100*t60 - t105*t60 - t107*t64 - t166*t99 + t172 - t28*(t128 + t184));
   double t186 = t116*t179 + t35*(-t117*t166 - t118*t60 - t122*t60 - t123*t64 - t28*(t170 + t84));
   double t187 = t133*t179 + t35*(-t130*t70 - t134*t60 + t135*t60 + t136*t64 + t156);
   double t188 = t143*t179 + t35*(-t144*t60 + t146*t60 + t147*t64 - t166*t5 - t177);
   double t189 = t154*t179 + t35*(-t155*t60 + t157*t60 + t158*t64 - t166*t8);
   double t190 = -t76 - t80;
   double t191 = -2*pow(y0, 2);
   double t192 = -2*pow(z0, 2);
   double t193 = t17*t26*t27*x0;
   double t194 = t190*t98 + t35*(-t100*t75 - t105*t75 - t107*t79 - t193*y0 - t86*t99);
   double t195 = t116*t190 + t35*(-t117*t86 - t118*t75 - t122*t75 - t123*t79 - t193*z0);
   double t196 = t133*t190 + t35*(-t134*t75 + t135*t75 + t136*t79 - t2*t86);
   double t197 = t101*z0;
   double t198 = t143*t190 + t35*(-t144*t75 + t146*t75 + t147*t79 - t197 - t5*t86);
   double t199 = t101*y0;
   double t200 = t154*t190 + t35*(-t155*t75 + t157*t75 + t158*t79 + t199 - t8*t86);
   double t201 = -t94 - t97;
   double t202 = -2*pow(x0, 2);
   double t203 = t17*t27*t9*t96;
   double t204 = t114*t17*t27*t9;
   double t205 = t116*t201 + t35*(-t117*t203 - t118*t93 - t123*t96 - t165*y0*z0 - t204*t99);
   double t206 = t133*t201 + t35*(-t105*t130 - t134*t93 + t135*t93 + t136*t96 + t197);
   double t207 = t143*t201 + t35*(-t144*t93 + t146*t93 + t147*t96 - t203*t5);
   double t208 = t101*x0;
   double t209 = t154*t201 + t35*(-t155*t93 + t157*t93 + t158*t96 - t203*t8 - t208);
   double t210 = -t113 - t115;
   double t211 = t133*t210 + t35*(t111*t135 - t112*t2 + t114*t136 - t122*t130 - t199);
   double t212 = t143*t210 + t35*(t111*t146 - t112*t5 + t114*t147 - t204*t5 + t208);
   double t213 = t154*t210 + t35*(-t111*t155 + t111*t157 + t114*t158 - t204*t8);
   double t214 = 1.0/(t16*t16*t16);
   double t215 = 4*t20*t214*t26;
   double t216 = grad*t32*(t127 - t131);
   double t217 = 4*t20*t214*t26*x2;
   double t218 = t142*t216 + t35*(-t130*t144 + t130*t146 + t135*t140 - t217*y2);
   double t219 = t153*t216 + t35*(-t130*t155 + t130*t157 + t135*t151 - t217*z2);
   double t220 = grad*t32*(t138 - t141);
   double t221 = t153*t220 + t35*(-t140*t155 + t140*t157 + t146*t151 - t215*y2*z2);
   // Fill up Hessian
   hessian[(3*i0+0)*ndelta + 3*i0+0] +=  t34*(t22 + t29) + t35*(-t12*t39 - t21*t36 - t28*(t37 + t38) - t41*(-4*t23 - 4*t24));
   hessian[(3*i0+0)*ndelta + 3*i0+1] +=  t57;
   hessian[(3*i0+0)*ndelta + 3*i0+2] +=  t72;
   hessian[(3*i0+1)*ndelta + 3*i0+0] +=  t57;
   hessian[(3*i0+1)*ndelta + 3*i0+1] +=  t161*t51 + t35*(-t163*t44 - t164*t56 - t28*(t162 + t38) - t44*t53);
   hessian[(3*i0+1)*ndelta + 3*i0+2] +=  t167;
   hessian[(3*i0+2)*ndelta + 3*i0+0] +=  t72;
   hessian[(3*i0+2)*ndelta + 3*i0+1] +=  t167;
   hessian[(3*i0+2)*ndelta + 3*i0+2] +=  t180*t66 + t35*(-t166*t60 - t181*t71 - t28*(t162 + t37) - t60*t69);
   hessian[(3*i0+0)*ndelta + 3*i1+0] +=  t90;
   hessian[(3*i0+0)*ndelta + 3*i1+1] +=  t108;
   hessian[(3*i0+0)*ndelta + 3*i1+2] +=  t124;
   hessian[(3*i0+1)*ndelta + 3*i1+0] +=  t169;
   hessian[(3*i0+1)*ndelta + 3*i1+1] +=  t171;
   hessian[(3*i0+1)*ndelta + 3*i1+2] +=  t174;
   hessian[(3*i0+2)*ndelta + 3*i1+0] +=  t183;
   hessian[(3*i0+2)*ndelta + 3*i1+1] +=  t185;
   hessian[(3*i0+2)*ndelta + 3*i1+2] +=  t186;
   hessian[(3*i0+0)*ndelta + 3*i2+0] +=  t137;
   hessian[(3*i0+0)*ndelta + 3*i2+1] +=  t148;
   hessian[(3*i0+0)*ndelta + 3*i2+2] +=  t159;
   hessian[(3*i0+1)*ndelta + 3*i2+0] +=  t175;
   hessian[(3*i0+1)*ndelta + 3*i2+1] +=  t176;
   hessian[(3*i0+1)*ndelta + 3*i2+2] +=  t178;
   hessian[(3*i0+2)*ndelta + 3*i2+0] +=  t187;
   hessian[(3*i0+2)*ndelta + 3*i2+1] +=  t188;
   hessian[(3*i0+2)*ndelta + 3*i2+2] +=  t189;
   hessian[(3*i1+0)*ndelta + 3*i0+0] +=  t90;
   hessian[(3*i1+0)*ndelta + 3*i0+1] +=  t169;
   hessian[(3*i1+0)*ndelta + 3*i0+2] +=  t183;
   hessian[(3*i1+1)*ndelta + 3*i0+0] +=  t108;
   hessian[(3*i1+1)*ndelta + 3*i0+1] +=  t171;
   hessian[(3*i1+1)*ndelta + 3*i0+2] +=  t185;
   hessian[(3*i1+2)*ndelta + 3*i0+0] +=  t124;
   hessian[(3*i1+2)*ndelta + 3*i0+1] +=  t174;
   hessian[(3*i1+2)*ndelta + 3*i0+2] +=  t186;
   hessian[(3*i1+0)*ndelta + 3*i1+0] +=  grad*t190*t32*t81 + t35*(-1.0L/2.0L*t17*t26*t40*t79*t89 - t28*(t191 + t192) - t75*t83 - t75*t86);
   hessian[(3*i1+0)*ndelta + 3*i1+1] +=  t194;
   hessian[(3*i1+0)*ndelta + 3*i1+2] +=  t195;
   hessian[(3*i1+1)*ndelta + 3*i1+0] +=  t194;
   hessian[(3*i1+1)*ndelta + 3*i1+1] +=  t201*t98 + t35*(-t100*t93 - t107*t96 - t203*t93 - t28*(t192 + t202));
   hessian[(3*i1+1)*ndelta + 3*i1+2] +=  t205;
   hessian[(3*i1+2)*ndelta + 3*i1+0] +=  t195;
   hessian[(3*i1+2)*ndelta + 3*i1+1] +=  t205;
   hessian[(3*i1+2)*ndelta + 3*i1+2] +=  t116*t210 + t35*(-t111*t204 - t112*t117 - t114*t123 - t28*(t191 + t202));
   hessian[(3*i1+0)*ndelta + 3*i2+0] +=  t196;
   hessian[(3*i1+0)*ndelta + 3*i2+1] +=  t198;
   hessian[(3*i1+0)*ndelta + 3*i2+2] +=  t200;
   hessian[(3*i1+1)*ndelta + 3*i2+0] +=  t206;
   hessian[(3*i1+1)*ndelta + 3*i2+1] +=  t207;
   hessian[(3*i1+1)*ndelta + 3*i2+2] +=  t209;
   hessian[(3*i1+2)*ndelta + 3*i2+0] +=  t211;
   hessian[(3*i1+2)*ndelta + 3*i2+1] +=  t212;
   hessian[(3*i1+2)*ndelta + 3*i2+2] +=  t213;
   hessian[(3*i2+0)*ndelta + 3*i0+0] +=  t137;
   hessian[(3*i2+0)*ndelta + 3*i0+1] +=  t175;
   hessian[(3*i2+0)*ndelta + 3*i0+2] +=  t187;
   hessian[(3*i2+1)*ndelta + 3*i0+0] +=  t148;
   hessian[(3*i2+1)*ndelta + 3*i0+1] +=  t176;
   hessian[(3*i2+1)*ndelta + 3*i0+2] +=  t188;
   hessian[(3*i2+2)*ndelta + 3*i0+0] +=  t159;
   hessian[(3*i2+2)*ndelta + 3*i0+1] +=  t178;
   hessian[(3*i2+2)*ndelta + 3*i0+2] +=  t189;
   hessian[(3*i2+0)*ndelta + 3*i1+0] +=  t196;
   hessian[(3*i2+0)*ndelta + 3*i1+1] +=  t206;
   hessian[(3*i2+0)*ndelta + 3*i1+2] +=  t211;
   hessian[(3*i2+1)*ndelta + 3*i1+0] +=  t198;
   hessian[(3*i2+1)*ndelta + 3*i1+1] +=  t207;
   hessian[(3*i2+1)*ndelta + 3*i1+2] +=  t212;
   hessian[(3*i2+2)*ndelta + 3*i1+0] +=  t200;
   hessian[(3*i2+2)*ndelta + 3*i1+1] +=  t209;
   hessian[(3*i2+2)*ndelta + 3*i1+2] +=  t213;
   hessian[(3*i2+0)*ndelta + 3*i2+0] +=  t132*t216 + t35*(t126 - t13*t215 - t130*t134 + t135*(t173 + t184));
   hessian[(3*i2+0)*ndelta + 3*i2+1] +=  t218;
   hessian[(3*i2+0)*ndelta + 3*i2+2] +=  t219;
   hessian[(3*i2+1)*ndelta + 3*i2+0] +=  t218;
   hessian[(3*i2+1)*ndelta + 3*i2+1] +=  t142*t220 + t35*(t126 - t14*t215 - t140*t144 + t146*(t120 + t182));
   hessian[(3*i2+1)*ndelta + 3*i2+2] +=  t221;
   hessian[(3*i2+2)*ndelta + 3*i2+0] +=  t219;
   hessian[(3*i2+2)*ndelta + 3*i2+1] +=  t221;
   hessian[(3*i2+2)*ndelta + 3*i2+2] +=  t154*(t149 - t152) + t35*(t126 - t15*t215 - t151*t155 + t157*(t103 + t168));
}


void hessian_oop_cos(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double* hessian, long ndelta) {
  dlist_row_type *delta0, *delta1, *delta2;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  hessian_oop_cos_low(ic, (*ic).i0, (*ic).i1, (*ic).i2, delta0, delta1, delta2, value, grad, hessian, ndelta);
}

void hessian_oop_meancos(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double* hessian, long ndelta) {
  double tmp;
  dlist_row_type *delta0, *delta1, *delta2;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  tmp = forward_oop_cos_low((double*)delta0, (double*)delta1, (double*)delta2);
  hessian_oop_cos_low(ic, (*ic).i0, (*ic).i1, (*ic).i2, delta0, delta1, delta2, tmp, grad/3.0, hessian, ndelta);
  tmp = forward_oop_cos_low((double*)delta2, (double*)delta0, (double*)delta1);
  hessian_oop_cos_low(ic, (*ic).i2, (*ic).i0, (*ic).i1, delta2, delta0, delta1, tmp, grad/3.0, hessian, ndelta);
  tmp = forward_oop_cos_low((double*)delta1, (double*)delta2, (double*)delta0);
  hessian_oop_cos_low(ic, (*ic).i1, (*ic).i2, (*ic).i0, delta1, delta2, delta0, tmp, grad/3.0, hessian, ndelta);
}

void hessian_oop_angle_low(iclist_row_type* ic, long i0, long i1, long i2, dlist_row_type *deltas, dlist_row_type *delta0, dlist_row_type *delta1, dlist_row_type *delta2, double value, double grad, double* hessian, long ndelta) {
   double tmp = sin(value);
   if (tmp!=0.0) tmp = -grad/tmp;
   // First part
   hessian_oop_cos_low(ic, i0, i1, i2, delta0, delta1, delta2, cos(value), tmp, hessian, ndelta);
   // Second part
   double jacobian[9];
   //delta0 = deltas + (*ic).i0;
   //delta1 = deltas + (*ic).i1;
   //delta2 = deltas + (*ic).i2;
   (*delta0).gx = 0.0;
   (*delta0).gy = 0.0;
   (*delta0).gz = 0.0;
   (*delta1).gx = 0.0;
   (*delta1).gy = 0.0;
   (*delta1).gz = 0.0;
   (*delta2).gx = 0.0;
   (*delta2).gy = 0.0;
   (*delta2).gz = 0.0;
   back_oop_angle_low(delta0, delta1, delta2, value, 1.0);
   jacobian[0] = (*delta0).gx;
   jacobian[1] = (*delta0).gy;
   jacobian[2] = (*delta0).gz;
   jacobian[3] = (*delta1).gx;
   jacobian[4] = (*delta1).gy;
   jacobian[5] = (*delta1).gz;
   jacobian[6] = (*delta2).gx;
   jacobian[7] = (*delta2).gy;
   jacobian[8] = (*delta2).gz;
   long i,j;
   tmp *= -cos(value);
   for (i=0;i<3;i++){
     for (j=0;j<3;j++){
       hessian[(3*i0+i)*ndelta + 3*i0+j] -= tmp*jacobian[i]*jacobian[j];
       hessian[(3*i0+i)*ndelta + 3*i1+j] -= tmp*jacobian[i]*jacobian[j+3];
       hessian[(3*i0+i)*ndelta + 3*i2+j] -= tmp*jacobian[i]*jacobian[j+6];
       hessian[(3*i1+i)*ndelta + 3*i0+j] -= tmp*jacobian[i+3]*jacobian[j];
       hessian[(3*i1+i)*ndelta + 3*i1+j] -= tmp*jacobian[i+3]*jacobian[j+3];
       hessian[(3*i1+i)*ndelta + 3*i2+j] -= tmp*jacobian[i+3]*jacobian[j+6];
       hessian[(3*i2+i)*ndelta + 3*i0+j] -= tmp*jacobian[i+6]*jacobian[j];
       hessian[(3*i2+i)*ndelta + 3*i1+j] -= tmp*jacobian[i+6]*jacobian[j+3];
       hessian[(3*i2+i)*ndelta + 3*i2+j] -= tmp*jacobian[i+6]*jacobian[j+6];
     }
   }
}

void hessian_oop_angle(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double* hessian, long ndelta) {
  dlist_row_type *delta0, *delta1, *delta2;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  hessian_oop_angle_low(ic, (*ic).i0, (*ic).i1, (*ic).i2, deltas, delta0, delta1, delta2, value, grad, hessian, ndelta);
}

void hessian_oop_meanangle(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double* hessian, long ndelta) {
  double tmp;
  dlist_row_type *delta0, *delta1, *delta2;
  delta0 = deltas + (*ic).i0;
  delta1 = deltas + (*ic).i1;
  delta2 = deltas + (*ic).i2;
  tmp = forward_oop_angle_low((double*)delta0, (double*)delta1, (double*)delta2);
  hessian_oop_angle_low(ic, (*ic).i0, (*ic).i1, (*ic).i2, deltas, delta0, delta1, delta2, tmp, grad/3.0, hessian, ndelta);
  tmp = forward_oop_angle_low((double*)delta2, (double*)delta0, (double*)delta1);
  hessian_oop_angle_low(ic, (*ic).i2, (*ic).i0, (*ic).i1, deltas, delta2, delta0, delta1, tmp, grad/3.0, hessian, ndelta);
  tmp = forward_oop_angle_low((double*)delta1, (double*)delta2, (double*)delta0);
  hessian_oop_angle_low(ic, (*ic).i1, (*ic).i2, (*ic).i0, deltas, delta1, delta2, delta0, tmp, grad/3.0, hessian, ndelta);
}

void hessian_oop_distance(iclist_row_type* ic, dlist_row_type* deltas, double value, double grad, double* hessian, long ndelta) {
   // Initialization
   dlist_row_type *delta0;
   delta0 = deltas + (*ic).i0;
   double x0 = (*delta0).dx;
   double y0 = (*delta0).dy;
   double z0 = (*delta0).dz;
   dlist_row_type *delta1;
   delta1 = deltas + (*ic).i1;
   double x1 = (*delta1).dx;
   double y1 = (*delta1).dy;
   double z1 = (*delta1).dz;
   dlist_row_type *delta2;
   delta2 = deltas + (*ic).i2;
   double x2 = (*delta2).dx;
   double y2 = (*delta2).dy;
   double z2 = (*delta2).dz;
   // Compute auxiliary variables
   double t0 = x0*y1;
   double t1 = x1*y0;
   double t2 = -t1;
   double t3 = t0 + t2;
   double t4 = x0*z1;
   double t5 = x1*z0;
   double t6 = -t5;
   double t7 = t4 + t6;
   double t8 = y0*z1;
   double t9 = y1*z0;
   double t10 = -t9;
   double t11 = t10 + t8;
   double t12 = t11*t11 + t3*t3 + t7*t7;
   double t13 = grad/sqrt(pow(t12, 5));
   double t14 = t3*y1 + t7*z1;
   double t15 = 3*t11*x2 + 3*t3*z2 - 3*t7*y2;
   double t16 = y1*z2 - y2*z1;
   double t17 = y1*y1;
   double t18 = z1*z1;
   double t19 = t11*x2 + t3*z2 - t7*y2;
   double t20 = -t11*z1 + t3*x1;
   double t21 = t14*t15;
   double t22 = x1*z2 - x2*z1;
   double t23 = t19*x1;
   double t24 = t13*(t12*(t14*t22 + t16*t20 + t23*y1) - t20*t21);
   double t25 = t11*y1 + t7*x1;
   double t26 = t15*t25;
   double t27 = x1*y2 - x2*y1;
   double t28 = t13*(t12*(-t14*t27 + t16*t25 + t23*z1) - t14*t26);
   double t29 = t3*y0 + t7*z0;
   double t30 = t15*t29;
   double t31 = y0*z2 - y2*z0;
   double t32 = y0*y1;
   double t33 = z0*z1;
   double t34 = t13*(t12*(t14*t31 + t16*t29 + t19*(t32 + t33)) - t14*t30);
   double t35 = t12*t12;
   double t36 = t35*z2;
   double t37 = -t11*z0 + t3*x0;
   double t38 = x0*z2 - x2*z0;
   double t39 = t13*(-t12*(t14*t38 + t16*t37 + t19*(2*t0 + t2)) + t21*t37 + t36);
   double t40 = t35*y2;
   double t41 = t11*y0 + t7*x0;
   double t42 = t15*t41;
   double t43 = x0*y2 - x2*y0;
   double t44 = t13*(t12*(t14*t43 - t16*t41 - t19*(2*t4 + t6)) + t14*t42 - t40);
   double t45 = 1.0/(sqrt(t12*t12*t12));
   double t46 = grad*t11*t45;
   double t47 = -t14*t46;
   double t48 = grad*t45;
   double t49 = t12*z1;
   double t50 = t48*(t14*t7 - t49);
   double t51 = t12*y1;
   double t52 = t48*(-t14*t3 + t51);
   double t53 = x1*x1;
   double t54 = t13*(t12*(t19*y1*z1 + t20*t27 - t22*t25) + t20*t26);
   double t55 = t13*(t12*(t19*(t0 - 2*t1) - t20*t31 - t22*t29) + t20*t30 - t36);
   double t56 = x0*x1;
   double t57 = t13*(t12*(t19*(t33 + t56) + t20*t38 + t22*t37) - t15*t20*t37);
   double t58 = t35*x2;
   double t59 = t13*(t12*(-t19*(t10 + 2*t8) - t20*t43 + t22*t41) - t20*t42 + t58);
   double t60 = t48*(t11*t20 + t49);
   double t61 = grad*t45*t7;
   double t62 = -t20*t61;
   double t63 = t12*x1;
   double t64 = t48*(t20*t3 - t63);
   double t65 = t13*(t12*(t19*(t4 - 2*t5) - t25*t31 + t27*t29) + t26*t29 + t40);
   double t66 = t13*(t12*(t19*(t8 - 2*t9) + t25*t38 - t27*t37) - t26*t37 - t58);
   double t67 = t13*(t12*(t19*(t32 + t56) - t25*t43 - t27*t41) - t25*t42);
   double t68 = t48*(t11*t25 - t51);
   double t69 = t48*(-t25*t7 + t63);
   double t70 = grad*t3*t45;
   double t71 = t25*t70;
   double t72 = y0*y0;
   double t73 = z0*z0;
   double t74 = t19*x0;
   double t75 = t13*(t12*(t29*t38 + t31*t37 + t74*y0) - t30*t37);
   double t76 = t13*(t12*(-t29*t43 + t31*t41 + t74*z0) - t29*t42);
   double t77 = t29*t46;
   double t78 = t12*z0;
   double t79 = t48*(-t29*t7 + t78);
   double t80 = t12*y0;
   double t81 = t48*(t29*t3 - t80);
   double t82 = x0*x0;
   double t83 = t13*(t12*(t19*y0*z0 + t37*t43 - t38*t41) + t37*t42);
   double t84 = -t48*(t11*t37 + t78);
   double t85 = t37*t61;
   double t86 = t12*x0;
   double t87 = t48*(-t3*t37 + t86);
   double t88 = t48*(-t11*t41 + t80);
   double t89 = t48*(t41*t7 - t86);
   double t90 = -t41*t70;
   // Fill up Hessian
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i0+0] +=  t13*(-t12*(2*t14*t16 + t19*(t17 + t18)) + pow(t14, 2)*t15);
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i0+1] +=  t24;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i0+2] +=  t28;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i0+0] +=  t24;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i0+1] +=  t13*(-t12*(t19*(t18 + t53) + 2*t20*t22) + t15*pow(t20, 2));
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i0+2] +=  t54;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i0+0] +=  t28;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i0+1] +=  t54;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i0+2] +=  t13*(t12*(-t19*(t17 + t53) + 2*t25*t27) + t15*pow(t25, 2));
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i1+0] +=  t34;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i1+1] +=  t39;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i1+2] +=  t44;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i1+0] +=  t55;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i1+1] +=  t57;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i1+2] +=  t59;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i1+0] +=  t65;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i1+1] +=  t66;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i1+2] +=  t67;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i2+0] +=  t47;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i2+1] +=  t50;
   hessian[(3*(*ic).i0+0)*ndelta + 3*(*ic).i2+2] +=  t52;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i2+0] +=  t60;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i2+1] +=  t62;
   hessian[(3*(*ic).i0+1)*ndelta + 3*(*ic).i2+2] +=  t64;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i2+0] +=  t68;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i2+1] +=  t69;
   hessian[(3*(*ic).i0+2)*ndelta + 3*(*ic).i2+2] +=  t71;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i0+0] +=  t34;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i0+1] +=  t55;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i0+2] +=  t65;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i0+0] +=  t39;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i0+1] +=  t57;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i0+2] +=  t66;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i0+0] +=  t44;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i0+1] +=  t59;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i0+2] +=  t67;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i1+0] +=  t13*(-t12*(t19*(t72 + t73) + 2*t29*t31) + t15*pow(t29, 2));
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i1+1] +=  t75;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i1+2] +=  t76;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i1+0] +=  t75;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i1+1] +=  t13*(-t12*(t19*(t73 + t82) + 2*t37*t38) + t15*pow(t37, 2));
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i1+2] +=  t83;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i1+0] +=  t76;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i1+1] +=  t83;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i1+2] +=  t13*(t12*(-t19*(t72 + t82) + 2*t41*t43) + t15*pow(t41, 2));
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i2+0] +=  t77;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i2+1] +=  t79;
   hessian[(3*(*ic).i1+0)*ndelta + 3*(*ic).i2+2] +=  t81;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i2+0] +=  t84;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i2+1] +=  t85;
   hessian[(3*(*ic).i1+1)*ndelta + 3*(*ic).i2+2] +=  t87;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i2+0] +=  t88;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i2+1] +=  t89;
   hessian[(3*(*ic).i1+2)*ndelta + 3*(*ic).i2+2] +=  t90;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i0+0] +=  t47;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i0+1] +=  t60;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i0+2] +=  t68;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i0+0] +=  t50;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i0+1] +=  t62;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i0+2] +=  t69;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i0+0] +=  t52;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i0+1] +=  t64;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i0+2] +=  t71;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i1+0] +=  t77;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i1+1] +=  t84;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i1+2] +=  t88;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i1+0] +=  t79;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i1+1] +=  t85;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i1+2] +=  t89;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i1+0] +=  t81;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i1+1] +=  t87;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i1+2] +=  t90;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i2+0] +=  0;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i2+1] +=  0;
   hessian[(3*(*ic).i2+0)*ndelta + 3*(*ic).i2+2] +=  0;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i2+0] +=  0;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i2+1] +=  0;
   hessian[(3*(*ic).i2+1)*ndelta + 3*(*ic).i2+2] +=  0;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i2+0] +=  0;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i2+1] +=  0;
   hessian[(3*(*ic).i2+2)*ndelta + 3*(*ic).i2+2] +=  0;
}


ic_hessian_type ic_hessian_fns[11] = {
  hessian_bond, hessian_bend_cos, hessian_bend_angle, hessian_dihed_cos, hessian_dihed_angle, hessian_bond,
  hessian_oop_cos, hessian_oop_meancos, hessian_oop_angle, hessian_oop_meanangle, hessian_oop_distance
};

void iclist_hessian(dlist_row_type* deltas, iclist_row_type* ictab, long nic, long ndelta, double* hessian) {
  long i;
  for (i=0; i<nic; i++) {
    ic_hessian_fns[ictab[i].kind](ictab + i, deltas, ictab[i].value, ictab[i].grad, hessian, ndelta);
  }
}
