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


#include "switchon.h"
#include "constants.h"
#include <math.h>
#include <stdlib.h>


switchon_type* switchon_new(void) {
  switchon_type* result;
  result = (switchon_type*)malloc(sizeof(switchon_type));
  if (result != NULL) {
    (*result).switchon_data = NULL;
    (*result).switchon_fn = NULL;
  }
  return result;
}

void switchon_free(switchon_type *switchon) {
  free(switchon);
}

int switchon_ready(switchon_type *switchon) {
  return (*switchon).switchon_data != NULL && (*switchon).switchon_fn != NULL;
}

void switchon_data_erf_init(switchon_type *switchon, double *radii) {
  switchon_data_erf_type *switchon_data;
  switchon_data = (switchon_data_erf_type*)malloc(sizeof(switchon_data_erf_type));
  (*switchon).switchon_data = switchon_data;
  if (switchon_data != NULL) {
    (*switchon).switchon_fn = switchon_fn_erf;
    (*switchon_data).radii = radii;
  }
}

double switchon_fn_erf(void *switchon_data, long center_index, long other_index, double d, double *g, double *gg) {
  double r_ab, y;
  r_ab = sqrt( (*(switchon_data_erf_type*)switchon_data).radii[center_index] * (*(switchon_data_erf_type*)switchon_data).radii[center_index] +
               (*(switchon_data_erf_type*)switchon_data).radii[other_index]  * (*(switchon_data_erf_type*)switchon_data).radii[other_index] );
  y = d/r_ab;
  if (r_ab > 0) {
    if (g != NULL) *g = M_TWO_DIV_SQRT_PI*exp(-y*y)/r_ab;
    if (gg != NULL) *gg = -2.0*M_TWO_DIV_SQRT_PI*y/(r_ab*r_ab)*exp(-y*y);
  }
  return erf(y);
}
