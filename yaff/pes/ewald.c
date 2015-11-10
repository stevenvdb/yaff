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
#include <stdlib.h>
#include "constants.h"
#include "ewald.h"
#include "cell.h"
#include <stdio.h>

double compute_ewald_reci(double *pos, long natom, double *charges,
                          cell_type* cell, double alpha, long *gmax, double
                          gcut, double dielectric, double *gpos, double *work,
                          double* vtens, double* hess) {
  long g0, g1, g2, i, j;
  double energy, k[3], ksq, cosfac, sinfac, x, c, s, h, fac1, fac2, dielectric_factor;
  double kvecs[9];
  for (i=0; i<9; i++) {
    kvecs[i] = M_TWO_PI*(*cell).gvecs[i];
  }
  energy = 0.0;
  fac1 = M_FOUR_PI/(*cell).volume;
  fac2 = 0.25/alpha/alpha;
  gcut *= M_TWO_PI;
  gcut *= gcut;
  for (g0=-gmax[0]; g0 <= gmax[0]; g0++) {
    for (g1=-gmax[1]; g1 <= gmax[1]; g1++) {
      for (g2=0; g2 <= gmax[2]; g2++) {
        if (g2==0) {
          if (g1<0) continue;
          if ((g1==0)&&(g0<=0)) continue;
        }
        k[0] = (g0*kvecs[0] + g1*kvecs[3] + g2*kvecs[6]);
        k[1] = (g0*kvecs[1] + g1*kvecs[4] + g2*kvecs[7]);
        k[2] = (g0*kvecs[2] + g1*kvecs[5] + g2*kvecs[8]);
        ksq = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
        if (ksq > gcut) continue;
        cosfac = 0.0;
        sinfac = 0.0;
        for (i=0; i<natom; i++) {
          x = k[0]*pos[3*i] + k[1]*pos[3*i+1] + k[2]*pos[3*i+2];
          c = charges[i]*cos(x);
          s = charges[i]*sin(x);
          cosfac += c;
          sinfac += s;
          if (gpos != NULL) {
            work[2*i] = c;
            work[2*i+1] = -s;
          }
        }
        c = fac1*exp(-ksq*fac2)/ksq;
        s = (cosfac*cosfac+sinfac*sinfac);
        energy += c*s;
        if (gpos != NULL) {
          x = 2.0*c;
          cosfac *= x;
          sinfac *= x;
          for (i=0; i<natom; i++) {
            x = cosfac*work[2*i+1] + sinfac*work[2*i];
            gpos[3*i] += k[0]*x;
            gpos[3*i+1] += k[1]*x;
            gpos[3*i+2] += k[2]*x;
          }
        }
        if (vtens != NULL) {
          c *= 2.0*(1.0/ksq+fac2)*s;
          vtens[0] += c*k[0]*k[0];
          vtens[4] += c*k[1]*k[1];
          vtens[8] += c*k[2]*k[2];
          x = c*k[1]*k[0];
          vtens[1] += x;
          vtens[3] += x;
          x = c*k[2]*k[0];
          vtens[2] += x;
          vtens[6] += x;
          x = c*k[2]*k[1];
          vtens[5] += x;
          vtens[7] += x;
        }
        if (hess != NULL) {
          /*
          Contributions to all hessian elements are proportional to:
            q_i*q_j*(sin(k.r_i)*sin(k.r_j) + cos(k.r_i)*cos(k.r_j))
            - \delta_{i,j}*( \sum_n q_n*cos(k.r_n)*q_i*cos(k.r_i) + \sum_n q_n*sin(k.r_n)*q_i*sin(k.r_i) )
          Such terms have to be multiplied by a prefactor and by the appropriate components of k
          */
          c = 2.0*fac1*exp(-ksq*fac2)/ksq;
          for (i=0; i<natom; i++) {
            for (j=i+1; j<natom; j++) {
              x = (work[2*i+1]*work[2*j+1]+work[2*i]*work[2*j])*c;
              s = x*k[0];
              h = s*k[0];
              hess[9*natom*i+3*j            ] += h;
              hess[9*natom*j+3*i            ] += h;
              h = s*k[1];
              hess[9*natom*i+3*j+1          ] += h;
              hess[9*natom*j+3*i+1          ] += h;
              hess[9*natom*i+3*j+    3*natom] += h;
              hess[9*natom*j+3*i+    3*natom] += h;
              h = s*k[2];
              hess[9*natom*i+3*j+2          ] += h;
              hess[9*natom*j+3*i+2          ] += h;
              hess[9*natom*i+3*j+    6*natom] += h;
              hess[9*natom*j+3*i+    6*natom] += h;
              s = x*k[1];
              h = s*k[1];
              hess[9*natom*i+3*j+1+  3*natom] += h;
              hess[9*natom*j+3*i+1+  3*natom] += h;
              h = s*k[2];
              hess[9*natom*i+3*j+2+  3*natom] += h;
              hess[9*natom*i+3*j+1+  6*natom] += h;
              hess[9*natom*j+3*i+2+  3*natom] += h;
              hess[9*natom*j+3*i+1+  6*natom] += h;
              h = x*k[2]*k[2];
              hess[9*natom*i+3*j+2+  6*natom] += h;
              hess[9*natom*j+3*i+2+  6*natom] += h;
            }
            x = (work[2*i+1]*work[2*i+1]+work[2*i]*work[2*i])*c - (cosfac*work[2*i]-sinfac*work[2*i+1]);
            hess[9*natom*i+3*i+0] += x*k[0]*k[0];
            hess[9*natom*i+3*i+1] += x*k[0]*k[1];
            hess[9*natom*i+3*i+2] += x*k[0]*k[2];
            hess[9*natom*i+3*i+0+3*natom] += x*k[1]*k[0];
            hess[9*natom*i+3*i+1+3*natom] += x*k[1]*k[1];
            hess[9*natom*i+3*i+2+3*natom] += x*k[1]*k[2];
            hess[9*natom*i+3*i+0+6*natom] += x*k[2]*k[0];
            hess[9*natom*i+3*i+1+6*natom] += x*k[2]*k[1];
            hess[9*natom*i+3*i+2+6*natom] += x*k[2]*k[2];
          }
        }
      }
    }
  }
  if (vtens != NULL) {
    vtens[0] -= energy;
    vtens[4] -= energy;
    vtens[8] -= energy;
  }
  //Corrections for dielectric constant
  dielectric_factor = 1.0/dielectric;
  if (gpos != NULL) {
    for (i=0; i<(3*natom); i++) {
      gpos[i] *= dielectric_factor;
    }
  }
  if (vtens != NULL) {
    for (i=0; i<9; i++) {
    vtens[i] *= dielectric_factor;
    }
  }
  if (hess != NULL) {
    for (i=0; i<(3*natom*3*natom); i++) {
      hess[i] *= dielectric_factor;
    }
  }
  energy *= dielectric_factor;
  return energy;
}

//TODO: lot of code overlap with original Ewald
//At the moment the idea is to make separate code for systems with monopoles and dipoles.
//If it turns out that adding zero dipoles does not increase computational cost, this separate
//code should become the main.
double compute_ewald_reci_dd(double *pos, long natom, double *charges, double *dipoles,
                          cell_type* cell, double alpha, long *gmax,
                          double gcut, double *gpos, double *work,
                          double* vtens) {
  long g0, g1, g2, i;
  double energy, k[3], ksq, cosfac_dd[3], sinfac_dd[3], x, c, s, fac1, fac2;
  double cosfac, sinfac;
  double kvecs[9];
  for (i=0; i<9; i++) {
    kvecs[i] = M_TWO_PI*(*cell).gvecs[i];
  }
  energy = 0.0;
  fac1 = M_FOUR_PI/(*cell).volume;
  fac2 = 0.25/alpha/alpha;
  gcut *= M_TWO_PI;
  gcut *= gcut;
  for (g0=-gmax[0]; g0 <= gmax[0]; g0++) {
    for (g1=-gmax[1]; g1 <= gmax[1]; g1++) {
      for (g2=0; g2 <= gmax[2]; g2++) {
        if (g2==0) {
          if (g1<0) continue;
          if ((g1==0)&&(g0<=0)) continue;
        }
        k[0] = (g0*kvecs[0] + g1*kvecs[3] + g2*kvecs[6]);
        k[1] = (g0*kvecs[1] + g1*kvecs[4] + g2*kvecs[7]);
        k[2] = (g0*kvecs[2] + g1*kvecs[5] + g2*kvecs[8]);
        ksq = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
        if (ksq > gcut) continue;
        cosfac_dd[0] = 0.0;
        cosfac_dd[1] = 0.0;
        cosfac_dd[2] = 0.0;
        sinfac_dd[0] = 0.0;
        sinfac_dd[1] = 0.0;
        sinfac_dd[2] = 0.0;
        cosfac = 0.0;
        sinfac = 0.0;
        for (i=0; i<natom; i++) {
          x = k[0]*pos[3*i] + k[1]*pos[3*i+1] + k[2]*pos[3*i+2];
          c = charges[i]*cos(x) + (k[0]*dipoles[3*i+0] + k[1]*dipoles[3*i+1] + k[2]*dipoles[3*i+2])*sin(x);
          s = charges[i]*sin(x) - (k[0]*dipoles[3*i+0] + k[1]*dipoles[3*i+1] + k[2]*dipoles[3*i+2])*cos(x);
          cosfac += c;
          sinfac += s;
          if (gpos != NULL) {
            work[2*i+0] = charges[i]*cos(x) + (k[0]*dipoles[3*i+0] + k[1]*dipoles[3*i+1] + k[2]*dipoles[3*i+2])*sin(x);
            work[2*i+1] =-charges[i]*sin(x) + (k[0]*dipoles[3*i+0] + k[1]*dipoles[3*i+1] + k[2]*dipoles[3*i+2])*cos(x);
          }
          if (vtens != NULL){
              cosfac_dd[0] +=-dipoles[3*i+0]*sin(x);
              cosfac_dd[1] +=-dipoles[3*i+1]*sin(x);
              cosfac_dd[2] +=-dipoles[3*i+2]*sin(x);
              sinfac_dd[0] += dipoles[3*i+0]*cos(x);
              sinfac_dd[1] += dipoles[3*i+1]*cos(x);
              sinfac_dd[2] += dipoles[3*i+2]*cos(x);
          }
        }
        c = fac1*exp(-ksq*fac2)/ksq;
        s = (cosfac*cosfac+sinfac*sinfac);
        energy += c*s;
        if (gpos != NULL) {
          x = 2.0*c;
          cosfac *= x;
          sinfac *= x;
          for (i=0; i<natom; i++) {
            x = cosfac*work[2*i+1] + sinfac*work[2*i];
            gpos[3*i+0] += k[0]*x;
            gpos[3*i+1] += k[1]*x;
            gpos[3*i+2] += k[2]*x;
          }
        }
        if (vtens != NULL) {
          c *= 2.0*(1.0/ksq+fac2)*s;
          vtens[0] += c*k[0]*k[0] + cosfac_dd[0]*k[0]*cosfac + sinfac_dd[0]*k[0]*sinfac;
          vtens[4] += c*k[1]*k[1] + cosfac_dd[1]*k[1]*cosfac + sinfac_dd[1]*k[1]*sinfac;
          vtens[8] += c*k[2]*k[2] + cosfac_dd[2]*k[2]*cosfac + sinfac_dd[2]*k[2]*sinfac;
          x = c*k[1]*k[0];
          vtens[1] += x + cosfac_dd[0]*k[1]*cosfac + sinfac_dd[0]*k[1]*sinfac;
          vtens[3] += x + cosfac_dd[1]*k[0]*cosfac + sinfac_dd[1]*k[0]*sinfac;
          x = c*k[2]*k[0];
          vtens[2] += x + cosfac_dd[0]*k[2]*cosfac + sinfac_dd[0]*k[2]*sinfac;
          vtens[6] += x + cosfac_dd[2]*k[0]*cosfac + sinfac_dd[2]*k[0]*sinfac;
          x = c*k[2]*k[1];
          vtens[5] += x + cosfac_dd[1]*k[2]*cosfac + sinfac_dd[1]*k[2]*sinfac;
          vtens[7] += x + cosfac_dd[2]*k[1]*cosfac + sinfac_dd[2]*k[1]*sinfac;
        }
      }
    }
  }
  if (vtens != NULL) {
    vtens[0] -= energy;
    vtens[4] -= energy;
    vtens[8] -= energy;
  }
  return energy;
}

double compute_ewald_corr(double *pos, double *charges,
                          cell_type *unitcell, double alpha,
                          scaling_row_type *stab, long nstab, double dielectric,
                          double *gpos, double *vtens, double *hess, long natom) {
  long i, center_index, other_index, hsize;
  double energy, delta[3], d, x, g, gg, pot, fac, dielectric_factor, rmin2, h, hh;
  energy = 0.0;
  g = 0.0;
  gg = 0.0;
  hsize = 3*natom;
  // Self-interaction correction (no gpos or vtens or hess contribution)
  x = alpha/M_SQRT_PI;
  for (i = 0; i < natom; i++) {
    energy -= x*charges[i]*charges[i];
  }
  // Scaling corrections
  for (i = 0; i < nstab; i++) {
    center_index = stab[i].a;
    other_index = stab[i].b;
    delta[0] = pos[3*other_index    ] - pos[3*center_index    ];
    delta[1] = pos[3*other_index + 1] - pos[3*center_index + 1];
    delta[2] = pos[3*other_index + 2] - pos[3*center_index + 2];
    cell_mic(delta, unitcell);
    d = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
    x = alpha*d;
    pot = erf(x)/d;
    fac = (1-stab[i].scale)*charges[other_index]*charges[center_index];
    if ((gpos != NULL) || (vtens != NULL)) {
      rmin2 = 1.0/(d*d);
      g = -fac*(M_TWO_DIV_SQRT_PI*alpha*exp(-x*x) - pot)*rmin2;
      if (hess != NULL) {
        gg = -fac*(pot-M_TWO_DIV_SQRT_PI*alpha*(1.0+x*x)*exp(-x*x))*2.0*rmin2*rmin2;
      }
    }
    if (gpos != NULL) {
      x = delta[0]*g;
      gpos[3*other_index  ] += x;
      gpos[3*center_index   ] -= x;
      x = delta[1]*g;
      gpos[3*other_index+1] += x;
      gpos[3*center_index +1] -= x;
      x = delta[2]*g;
      gpos[3*other_index+2] += x;
      gpos[3*center_index +2] -= x;
    }
    if (vtens != NULL) {
      vtens[0] += delta[0]*delta[0]*g;
      vtens[4] += delta[1]*delta[1]*g;
      vtens[8] += delta[2]*delta[2]*g;
      x = delta[1]*delta[0]*g;
      vtens[1] += x;
      vtens[3] += x;
      x = delta[2]*delta[0]*g;
      vtens[2] += x;
      vtens[6] += x;
      x = delta[2]*delta[1]*g;
      vtens[5] += x;
      vtens[7] += x;
    }
    if (hess != NULL) {
        h = gg*delta[0];
        // x-x entries
        hh = h*delta[0] + g*(1.0-delta[0]*delta[0]*rmin2);
        hess[3*center_index*hsize+3*center_index+    0] += hh;
        hess[3*other_index*hsize +3*other_index +    0] += hh;
        hess[3*center_index*hsize+3*other_index +    0] -= hh;
        hess[3*other_index*hsize +3*center_index+    0] -= hh;
        // x-y and y-x entries
        hh = h*delta[1] - g*delta[0]*delta[1]*rmin2;
        hess[3*center_index*hsize+3*center_index+    1] += hh;
        hess[3*center_index*hsize+3*center_index+hsize] += hh;
        hess[3*other_index*hsize +3*other_index +    1] += hh;
        hess[3*other_index*hsize +3*other_index +hsize] += hh;
        hess[3*center_index*hsize+3*other_index +    1] -= hh;
        hess[3*center_index*hsize+3*other_index +hsize] -= hh;
        hess[3*other_index*hsize +3*center_index+    1] -= hh;
        hess[3*other_index*hsize +3*center_index+hsize] -= hh;

        h = gg*delta[1];
        // y-y entries
        hh = h*delta[1] + g*(1.0-delta[1]*delta[1]*rmin2);
        hess[3*center_index*hsize+3*center_index+hsize+1] += hh;
        hess[3*other_index*hsize +3*other_index +hsize+1] += hh;
        hess[3*center_index*hsize+3*other_index +hsize+1] -= hh;
        hess[3*other_index*hsize +3*center_index+hsize+1] -= hh;
        // y-z and z-y entries
        hh = h*delta[2] - g*delta[1]*delta[2]*rmin2;
        hess[3*center_index*hsize+3*center_index+hsize+2] += hh;
        hess[3*center_index*hsize+3*center_index+2*hsize+1] += hh;
        hess[3*other_index*hsize +3*other_index +hsize+2] += hh;
        hess[3*other_index*hsize +3*other_index +2*hsize+1] += hh;
        hess[3*center_index*hsize+3*other_index +hsize+2] -= hh;
        hess[3*center_index*hsize+3*other_index +2*hsize+1] -= hh;
        hess[3*other_index*hsize +3*center_index+hsize+2] -= hh;
        hess[3*other_index*hsize +3*center_index+2*hsize+1] -= hh;

        h = gg*delta[2];
        // z-z entries
        hh = h*delta[2] + g*(1.0-delta[2]*delta[2]*rmin2);
        hess[3*center_index*hsize+3*center_index+2*hsize+2] += hh;
        hess[3*other_index*hsize +3*other_index +2*hsize+2] += hh;
        hess[3*center_index*hsize+3*other_index +2*hsize+2] -= hh;
        hess[3*other_index*hsize +3*center_index+2*hsize+2] -= hh;
        // x-z and x-z entries
        hh = h*delta[0] - g*delta[2]*delta[0]*rmin2;
        hess[3*center_index*hsize+3*center_index+    2] += hh;
        hess[3*center_index*hsize+3*center_index+2*hsize] += hh;
        hess[3*other_index*hsize +3*other_index +    2] += hh;
        hess[3*other_index*hsize +3*other_index +2*hsize] += hh;
        hess[3*center_index*hsize+3*other_index +    2] -= hh;
        hess[3*center_index*hsize+3*other_index +2*hsize] -= hh;
        hess[3*other_index*hsize +3*center_index+    2] -= hh;
        hess[3*other_index*hsize +3*center_index+2*hsize] -= hh;
    }
    energy -= fac*pot;
  }
  //Corrections for dielectric constant
  dielectric_factor = 1.0/dielectric;
  if (gpos != NULL) {
    for (i=0; i<(3*natom); i++) {
      gpos[i] *= dielectric_factor;
    }
  }
  if (vtens != NULL) {
    for (i=0; i<9; i++) {
    vtens[i] *= dielectric_factor;
    }
  }
  if (hess != NULL) {
    for (i=0; i<(3*natom*3*natom); i++) {
      hess[i] *= dielectric_factor;
    }
  }
  energy *= dielectric_factor;
  return energy;
}


double compute_ewald_corr_dd(double *pos, double *charges, double *dipoles,
                          cell_type *unitcell, double alpha,
                          scaling_row_type *stab, long nstab,
                          double *gpos, double *vtens, long natom) {
  long i, j, k;
  double energy, delta[3], d, x, g, g_cart[3];
  double fac, fac0, fac1, fac2, fac3, d_2;
  double mui_dot_delta, muj_dot_delta, mui_dot_muj;
  //double pot_cc, pot_cd, pot_dd;
  energy = 0.0;
  g = 0.0;
  fac1 = alpha/M_SQRT_PI;
  fac2 = fac1*alpha*alpha*2.0/3.0;
  // Self-interaction correction (no gpos or vtens contribution)
  for (i = 0; i < natom; i++) {
    //charges
    energy -= fac1*charges[i]*charges[i];
    //dipoles
    energy -= fac2*( dipoles[3*i+0]*dipoles[3*i+0] + dipoles[3*i+1]*dipoles[3*i+1] + dipoles[3*i+2]*dipoles[3*i+2] );
  }
  // Scaling corrections
  for (k = 0; k < nstab; k++) { // Loop over all pairs that need scaling
    i = stab[k].a;
    j = stab[k].b;
    delta[0] = pos[3*j+0] - pos[3*i+0];
    delta[1] = pos[3*j+1] - pos[3*i+1];
    delta[2] = pos[3*j+2] - pos[3*i+2];
    cell_mic(delta, unitcell);
    d = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
    //Some useful definitions
    d_2 = 1.0/(d*d);
    x = alpha*d;
    fac = (1-stab[i].scale);
    fac0 = erf(x)/d*fac;
    fac1 = (    fac0 - M_TWO_DIV_SQRT_PI*alpha*exp(-x*x)*fac)*d_2;
    fac2 = (3.0*fac1 - 2.0*M_TWO_DIV_SQRT_PI*alpha*alpha*alpha*exp(-x*x)*fac)*d_2;
    mui_dot_delta = dipoles[3*i+0]*delta[0] + dipoles[3*i+1]*delta[1] + dipoles[3*i+2]*delta[2];
    muj_dot_delta = dipoles[3*j+0]*delta[0] + dipoles[3*j+1]*delta[1] + dipoles[3*j+2]*delta[2];
    mui_dot_muj = dipoles[3*i+0]*dipoles[3*j+0] + dipoles[3*i+1]*dipoles[3*j+1] + dipoles[3*i+2]*dipoles[3*j+2];
    //CC interaction
    energy -= fac0*charges[j]*charges[i];
    //CD and DC interaction
    energy -= fac1*(charges[i]*muj_dot_delta - charges[j]*mui_dot_delta);
    //DD interaction
    energy -= (fac1*mui_dot_muj - fac2*mui_dot_delta*muj_dot_delta);
    if ((gpos != NULL) || (vtens != NULL)) {
      fac3 = (5.0*fac2 - 4.0*M_TWO_DIV_SQRT_PI*alpha*alpha*alpha*alpha*alpha*exp(-x*x)*fac)*d_2;
      //CC interaction
      g  = fac1*charges[j]*charges[i];
      //CD and DC interaction
      g += fac2*(charges[i]*muj_dot_delta - charges[j]*mui_dot_delta);
      //DD interaction
      g += fac2*mui_dot_muj - fac3*mui_dot_delta*muj_dot_delta;
      //CD and DC interaction
      g_cart[0] = fac1*(charges[j]*dipoles[3*i+0]-charges[i]*dipoles[3*j+0]);
      g_cart[1] = fac1*(charges[j]*dipoles[3*i+1]-charges[i]*dipoles[3*j+1]);
      g_cart[2] = fac1*(charges[j]*dipoles[3*i+2]-charges[i]*dipoles[3*j+2]);
      //DD interaction
      g_cart[0] += fac2*(dipoles[3*i+0]*muj_dot_delta + dipoles[3*j+0]*mui_dot_delta);
      g_cart[1] += fac2*(dipoles[3*i+1]*muj_dot_delta + dipoles[3*j+1]*mui_dot_delta);
      g_cart[2] += fac2*(dipoles[3*i+2]*muj_dot_delta + dipoles[3*j+2]*mui_dot_delta);
    }
    if (gpos != NULL) {
      x = delta[0]*g;
      gpos[3*j+0 ] += x + g_cart[0];
      gpos[3*i+0 ] -= x + g_cart[0];
      x = delta[1]*g;
      gpos[3*j+1 ] += x + g_cart[1];
      gpos[3*i+1 ] -= x + g_cart[1];
      x = delta[2]*g;
      gpos[3*j+2 ] += x + g_cart[2];
      gpos[3*i+2 ] -= x + g_cart[2];
    }
    if (vtens != NULL) {
      vtens[0] += delta[0]*(delta[0]*g+g_cart[0]);
      vtens[4] += delta[1]*(delta[1]*g+g_cart[1]);
      vtens[8] += delta[2]*(delta[2]*g+g_cart[2]);
      vtens[1] += delta[0]*(delta[1]*g+g_cart[1]);
      vtens[3] += delta[1]*(delta[0]*g+g_cart[0]);
      vtens[2] += delta[0]*(delta[2]*g+g_cart[2]);
      vtens[6] += delta[2]*(delta[0]*g+g_cart[0]);
      vtens[5] += delta[1]*(delta[2]*g+g_cart[2]);
      vtens[7] += delta[2]*(delta[1]*g+g_cart[1]);
    }
  }
  return energy;
}
