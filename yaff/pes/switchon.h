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


#ifndef YAFF_SWITCHON_H
#define YAFF_SWITCHON_H


typedef double (*switchon_fn_type)(void*, long, long, double, double*, double*);

typedef struct {
  void *switchon_data;
  switchon_fn_type switchon_fn;
} switchon_type;

switchon_type* switchon_new(void);
void switchon_free(switchon_type *switchon);
int switchon_ready(switchon_type *switchon);

typedef struct {
  double *radii;
} switchon_data_erf_type;

void switchon_data_erf_init(switchon_type *pair_pot, double *radii);
double switchon_fn_erf(void *switchon_data, long center_index, long other_index, double d, double *g, double *gg);


#endif
