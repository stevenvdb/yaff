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

#define DEBUG

#ifdef DEBUG
#include <stdio.h>
#endif

#include <math.h>
#include "nlist.h"
#include "cell.h"
#include "omp.h"


int nlist_build_low(double *pos, double rcut, long *rmax,
                    cell_type *unitcell, long *status,
                    neigh_row_type *neighs, long amax, long nneigh, long bmin, long bmax) {

  long a, b, row;
  long *r;
  int update_delta0, image, sign;
  double delta0[3], delta[3], d;

  // Compute square of the rcut distance
  rcut *= rcut;
  r = status;
  a = status[3];
  b = status[4];
  sign = status[5];

  update_delta0 = 1;
  image = (r[0] != 0) || (r[1] != 0) || (r[2] != 0);
  row = 0;

  while (row < nneigh) {
    if (a >= amax) {
      // Completely done.
      status[6] += row;
      return 1;
    }
    // Avoid adding pairs for which a > b and that match the minimum image
    // convention.
    if (update_delta0) {
      // Compute the relative vector.
      delta0[0] = pos[3*b  ] - pos[3*a  ];
      delta0[1] = pos[3*b+1] - pos[3*a+1];
      delta0[2] = pos[3*b+2] - pos[3*a+2];
      // Subtract the cell vectors as to make the relative vector as short
      // as possible. (This is the minimum image convention.)
      cell_mic(delta0, unitcell);
      // Done updating delta0.
      update_delta0 = 0;
    }
    // Only add self-interactions with atoms in periodic images.
    if ((b<a) || image) {
      // Construct delta by adding the appropriate cell vector to delta0
      delta[0] = sign*delta0[0];
      delta[1] = sign*delta0[1];
      delta[2] = sign*delta0[2];
      cell_add_vec(delta, unitcell, r);
      // Compute square of the distance and store the record if distance is below the rcut.
      // Square root is only computed if the record needs to be stored,
      // this way we avoid expensive square root in a lot of cases.
      d = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
      status[7]++;
      if (d < rcut) {
        if (sign > 0) {
          (*neighs).a = a;
          (*neighs).b = b;
        } else {
          (*neighs).a = b;
          (*neighs).b = a;
        }
        (*neighs).d = sqrt(d);
        (*neighs).dx = delta[0];
        (*neighs).dy = delta[1];
        (*neighs).dz = delta[2];
        (*neighs).r0 = r[0];
        (*neighs).r1 = r[1];
        (*neighs).r2 = r[2];
        neighs++;
        row++;
      }
    }
    // Increase the appropriate counters in the sextuple loop.
    if ((sign > 0) && (image) && (a!=b)) {
      // Change sign of the relative vector for non-self interactions with
      // periodic images.
      sign = -1;
    } else if (!nlist_inc_r(unitcell, r, rmax)) {
      sign = 1;
      update_delta0 = 1;
      image = 0;
      b++;
      if ((b > a) || (b > bmax)) {
      //if ((b > a)) {
        b = bmin;
        a++;
      }
    } else {
      image = 1;
      sign = 1;
    }
  }
  // Exit before the job is done. Keep track of the status. Work will be resumed
  // in a next call.
  status[0] = r[0];
  status[1] = r[1];
  status[2] = r[2];
  status[3] = a;
  status[4] = b;
  status[5] = sign;
  status[6] += row;
  return 0;
}


double decompose_domain(double *pos, cell_type *unitcell, long *bin_indexes, long *order, long natom, long *domains) {
  long i, j, index, bin, start, nneigh_total, starta, startb, stopa, stopb, nbins;
  double r_circum;
  long bins[natom];
  double small_rvecs[9], tmp[3], cart[3];
  // Determine number of domains along each cell vector
  //domains[0] = 1;
  //domains[1] = 1;
  //domains[2] = 1;
  // Initialize arrays that contain info about domains
  nbins = domains[0]*domains[1]*domains[2];
  long nbin[nbins], work[nbins], nbincum[nbins];
  for (i=0;i<nbins;i++){
    nbin[i] = 0;
    work[i] = 0;
    nbincum[i] = 0;
  }
  // Compute cell vectors of each domain
  for (i=0;i<9;i++) {
    small_rvecs[i] = (*unitcell).rvecs[i]/domains[i/3];
  }
  // Compute upper bound on the diameter of sphere circumscribing domain
  for (i=0;i<3;i++) tmp[i] = small_rvecs[0+i] + small_rvecs[3+i] + small_rvecs[6+i];
  r_circum = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
  // Put every atom in the correct domain
  for (i=0;i<natom;i++) {
    cart[0] = pos[3*i+0];
    cart[1] = pos[3*i+1];
    cart[2] = pos[3*i+2];
    cell_to_frac(unitcell, cart, tmp);
    bin_indexes[i]  = domains[1]*domains[2]*( (long) (domains[0]*(tmp[0] - floor(tmp[0]))));
    bin_indexes[i] += domains[2]*( (long) (domains[1]*(tmp[1] - floor(tmp[1]))));
    bin_indexes[i] += ( (long) (domains[2]*(tmp[2] - floor(tmp[2]))));
    nbin[bin_indexes[i]] += 1;
  }
  // Resort
  for (i=0;i<natom;i++) {
    bin = bin_indexes[i];
    index = work[bin];
    work[bin]++;
    for (j=0;j<bin;j++) index += nbin[j];
    order[index] = i;
  }
  return r_circum;
}


int nlist_domain_decomposition_omp(double *pos, double rcut, long *rmax,
                    cell_type *unitcell, long *status,
                    neigh_row_type *neighs, long natom, long nneigh, long nbin, long *binsizes, long *binsizes_cum,
                    long *domains, double r_circum, long nthreads) {
  long i,j, succes, startb, stopa, nneigh_total, ndist_total, nint, iint, failed, final_fail;
  omp_set_num_threads(nthreads);
  nneigh_total = 0;
  ndist_total = 0;

  double small_rvecs[9];
  long k, l, m, n;
  for (k=0;k<3;k++) {
     for (l=0;l<3;l++) {
        small_rvecs[3*k+l] = (*unitcell).rvecs[3*k+l]/domains[k];
        //printf("SMALL RVECS %d %6.2f\n",3*k+l,small_rvecs[3*k+l]);
     }
  }

  double r_centers;
  nint = (nbin*(nbin+1))/2;
  failed = 0;
  final_fail = 0;

  // Do the interactions inside different domains
  #pragma omp parallel private(i,j,k,l,m,succes, r_centers) firstprivate(nneigh_total, ndist_total, neighs, failed)
  {
  long nneigh_thread, thread_id;
  //nthreads = omp_get_num_threads();
  thread_id = omp_get_thread_num();
  nneigh_thread = nneigh/nthreads;
  neighs += thread_id*nneigh_thread;
  #pragma omp critical
  {
    //printf("This is thread %d out of %d, I will do at most %10d neighs, current pointer at %20d\n",thread_id, omp_get_num_threads(), nneigh_thread, neighs);
  }


  double delta[3];
  long mystatus[8];

  #pragma omp for
  for (iint=0;iint<nint;iint++) {
      i = (long) (sqrt(2*iint+0.25)-0.5);
      j = iint - (i*(i+1))/2;
      mystatus[0] = 0;
      mystatus[1] = 0;
      mystatus[2] = 0;
      mystatus[3] = binsizes_cum[i];
      mystatus[4] = binsizes_cum[j];
      mystatus[5] = 1;
      mystatus[6] = 0;
      mystatus[7] = 0;

        m = i%domains[2];
        l = (i/domains[2])%domains[1];
        k = (i/domains[2]/domains[1])%domains[0];
        delta[0] = k*small_rvecs[0] + l*small_rvecs[3] + m*small_rvecs[6];
        delta[1] = k*small_rvecs[1] + l*small_rvecs[4] + m*small_rvecs[7];
        delta[2] = k*small_rvecs[2] + l*small_rvecs[5] + m*small_rvecs[8];
        m = j%domains[2];
        l = (j/domains[2])%domains[1];
        k = (j/domains[2]/domains[1])%domains[0];
        delta[0] -= k*small_rvecs[0] + l*small_rvecs[3] + m*small_rvecs[6];
        delta[1] -= k*small_rvecs[1] + l*small_rvecs[4] + m*small_rvecs[7];
        delta[2] -= k*small_rvecs[2] + l*small_rvecs[5] + m*small_rvecs[8];
        cell_mic(delta, unitcell);
        r_centers = sqrt( delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2] );
      if (r_centers>rcut+r_circum || binsizes[i]==0 || binsizes[j]==0) {
          succes = 1;
      }
      else {
          // We only want interactions with one particle in domain "i" and the
          // particle in domain "j"
          stopa = binsizes_cum[i] + binsizes[i];
          startb = binsizes_cum[j];
          succes = nlist_build_low(pos, rcut, rmax,
                          unitcell, mystatus, neighs,
                          stopa, nneigh_thread-nneigh_total,
                          startb, startb+binsizes[j]-1);
      }
      //TODO Deal with not succes
      failed += (1-succes);
      //if (!succes) return 0;
      if (succes) {
      #ifdef DEBUG
      //printf("INTER %3d %3d by thread %d: %5d | %3d %3d | %d %d | %4.1f %4.1f %4.1f\n", i, j, thread_id, mystatus[6], stopa, startb, mystatus[3], mystatus[4], r_centers, r_circum, rcut);
      #endif
      // Update total number of neighbours found
      nneigh_total += mystatus[6];
      // Update total number of distance evaluations
      ndist_total += mystatus[7];
      // Increment pointer to neighs array with number of neighbours just found
      neighs += mystatus[6];}
  }
  #pragma omp critical
  {
    //printf("This is thread %d out of %d, I did %10d neighs and %10d distances, current pointer at %20d, failed %d\n",thread_id, nthreads, nneigh_total, ndist_total,  neighs, failed);
    status[3*thread_id] = thread_id*nneigh_thread;
    status[3*thread_id+1] = nneigh_total;
    status[3*thread_id+2] = ndist_total;
    final_fail += failed;
  }
  }
  if (final_fail) return 0;
  else return 1;
}



int nlist_domain_decomposition(double *pos, double rcut, long *rmax,
                    cell_type *unitcell, long *status,
                    neigh_row_type *neighs, long natom, long nneigh, long nbin, long *binsizes, long *binsizes_cum,
                    long *ijstart, long *domains, double r_circum) {
  long i,j, succes, startb, stopa, nneigh_total, ndist_total, istart, jstart;
  // Pick up where we left after previous build.
  istart = ijstart[0];
  jstart = ijstart[1];
  nneigh_total = 0;
  ndist_total = 0;
  // Do the interactions within same domain
  // jstart<0 indicates that previous build did not complete same-domain interactions
  if (jstart<0) {
      for (i=istart;i<nbin;i++){
        // We only want interactions within domain "i"
        stopa = binsizes_cum[i] + binsizes[i];
        startb = binsizes_cum[i];
        if (binsizes[i]==0) {
            succes = 1;
        } else {
            succes = nlist_build_low(pos, rcut, rmax,
                        unitcell, status, neighs,
                        stopa, nneigh-nneigh_total,
                        startb, startb+binsizes[i]-1);
        }
        #ifdef DEBUG
        printf("INTRA   %5d: %5d\n", i, status[6]);
        #endif
        // Update total number of neighbours found
        nneigh_total += status[6];
        // Update total number of distance evaluations
        ndist_total += status[7];
        // Increment pointer to neighs array with number of neighbours just found
        neighs += status[6];
        if (succes) {
            // This domain is completed! Reset status.
            status[0] = 0;
            status[1] = 0;
            status[2] = 0;
            status[5] = 1;
            status[6] = 0;
            status[7] = 0;
            if (nbin==1) {
                // No inter-domains to do!
                status[6] = nneigh_total;
                status[7] = ndist_total;
                return 1;
            }
            if (i==nbin-1) {
                // We are completely done with intra-domain interactions,
                // prepare to start with inter-domain interactions.
                istart = 1;
                jstart = 0;
                status[3] = binsizes_cum[istart];
                status[4] = binsizes_cum[jstart];
            }
            else {
              // We still have to do more intra-domain interactions,
              // prepare for the next domain
              status[3] = binsizes_cum[i+1];
              status[4] = binsizes_cum[i+1];
            }
        }
        else {
            // The neighs array was completely filled before all domains were completed,
            // save information to start next build where we are now.
            ijstart[0] = i;
            status[6] = nneigh_total;
            status[7] = ndist_total;
            return 0;
        }
      }
  }

  double small_rvecs[9];
  long k, l, m, n;
  for (k=0;k<3;k++) {
     for (l=0;l<3;l++) {
        small_rvecs[3*k+l] = (*unitcell).rvecs[3*k+l]/domains[k];
        //printf("SMALL RVECS %d %6.2f\n",3*k+l,small_rvecs[3*k+l]);
     }
  }
  double delta[3];
  double r_centers;

  // Do the interactions inside different domains, pick up where we left off
  // during previous build
  for (i=istart;i<nbin;i++){
    for (j=jstart;j<i;j++){
        m = i%domains[2];
        l = (i/domains[2])%domains[1];
        k = (i/domains[2]/domains[1])%domains[0];
        delta[0] = k*small_rvecs[0] + l*small_rvecs[3] + m*small_rvecs[6];
        delta[1] = k*small_rvecs[1] + l*small_rvecs[4] + m*small_rvecs[7];
        delta[2] = k*small_rvecs[2] + l*small_rvecs[5] + m*small_rvecs[8];
        m = j%domains[2];
        l = (j/domains[2])%domains[1];
        k = (j/domains[2]/domains[1])%domains[0];
        delta[0] -= k*small_rvecs[0] + l*small_rvecs[3] + m*small_rvecs[6];
        delta[1] -= k*small_rvecs[1] + l*small_rvecs[4] + m*small_rvecs[7];
        delta[2] -= k*small_rvecs[2] + l*small_rvecs[5] + m*small_rvecs[8];
        cell_mic(delta, unitcell);
        r_centers = sqrt( delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2] );
      if (r_centers>rcut+r_circum || binsizes[i]==0 || binsizes[j]==0) {
          succes = 1;
      }
      else {
          // We only want interactions with one particle in domain "i" and the
          // particle in domain "j"
          stopa = binsizes_cum[i] + binsizes[i];
          startb = binsizes_cum[j];
          succes = nlist_build_low(pos, rcut, rmax,
                          unitcell, status, neighs,
                          stopa, nneigh-nneigh_total,
                          startb, startb+binsizes[j]-1);
      }
      #ifdef DEBUG
      //printf("INTER %3d %3d by thread %d: %5d | %3d %3d | %d %d | %4.1f %4.1f %4.1f\n", i, j, thread_id, status[6], stopa, startb, status[3], status[4], r_centers, r_circum, rcut);
      printf("INTER %3d %3d: %5d | %3d %3d | %d %d | %4.1f %4.1f %4.1f\n", i, j, status[6], stopa, startb, status[3], status[4], r_centers, r_circum, rcut);
      #endif
      // Update total number of neighbours found
      nneigh_total += status[6];
      // Update total number of distance evaluations
      ndist_total += status[7];
      // Increment pointer to neighs array with number of neighbours just found
      neighs += status[6];
      if (succes) {
          // Interactions between domain "i" and "j" completed!
          // Reset status
          status[0] = 0;
          status[1] = 0;
          status[2] = 0;
          status[5] = 1;
          status[6] = 0;
          status[7] = 0;
          if ((i==nbin-1) && (j==i-1)) {
            // We are completely done! Store some final results in the status.
            status[6] = nneigh_total;
            status[7] = ndist_total;
            return 1;
          }
          else{
              // Prepare to do next inter-domain
              if (j==i-1) {
                jstart = 0;
                status[3] = binsizes_cum[i+1];
                status[4] = binsizes_cum[0];
              }
              else {
                status[3] = binsizes_cum[i];
                status[4] = binsizes_cum[j+1];;
              }
          }
      }
      else {
        // The neighs array was completely filled before all domains were completed,
        // save information to start next build where we are now.
        ijstart[0] = i;
        ijstart[1] = j;
        status[6] = nneigh_total;
        status[7] = ndist_total;
        return 0;
      }
    }
  }
}


int nlist_inc_r(cell_type *unitcell, long *r, long *rmax) {
  // increment the counters for the periodic images.
  // returns 1 when the counters were incremented successfully.
  // returns 0 and resets all the counters when the iteration over all cells
  // is complete.

  // Note: Only the central image and half of the neighboring images are
  // considered. This way, one can build neighborlists without any duplicates
  // pairs. The central cell comes first.
  if ((*unitcell).nvec > 0) {
    r[0]++;
    if (r[0] > rmax[0]) {
      r[0] = -rmax[0];
      if ((*unitcell).nvec > 1) {
        r[1]++;
        if (r[1] > rmax[1]) {
          r[1] = -rmax[1];
          if ((*unitcell).nvec > 2) {
            r[2]++;
            if (r[2] > rmax[2]) {
              r[2] = -rmax[2];
              goto endloop;
            }
          } else {
            goto endloop;
          }
        }
      } else {
        goto endloop;
      }
    }
  } else {
    return 0;
  }
  return 1;

endloop:
  // This point is only reached when all relevant neighboring cells are visited.
  // We set the counters back to the central cell
  if ((*unitcell).nvec > 0) r[0] = 0;
  if ((*unitcell).nvec > 1) r[1] = 0;
  if ((*unitcell).nvec > 2) r[2] = 0;
  return 0;
}


void nlist_recompute_low(double *pos, double *pos_old, cell_type* unitcell,
                         neigh_row_type *neighs, long nneigh) {
  long i, a, b;
  int update_delta0;
  long center[3];
  double delta0[3], delta[3], d;

  update_delta0 = 1;
  a = -1;
  b = -1;
  #pragma omp parallel for firstprivate(update_delta0, a, b) private(center, delta0, delta, d) ordered schedule(static)
  for (i=0; i<nneigh; i++) {
    if (neighs[i].a != a) {
      update_delta0 = 1;
    } else if (neighs[i].b != b) {
      update_delta0 = 1;
    }
    if (update_delta0) {
      a = neighs[i].a;
      b = neighs[i].b;
      // Compute the old relative vector.
      delta0[0] = pos_old[3*b  ] - pos_old[3*a  ];
      delta0[1] = pos_old[3*b+1] - pos_old[3*a+1];
      delta0[2] = pos_old[3*b+2] - pos_old[3*a+2];
      // Compute the cell vectors to be subtracted to bring the old to the MIC
      cell_to_center(delta0, unitcell, center);
      // Compute the new relative vector.
      delta0[0] = pos[3*b  ] - pos[3*a  ];
      delta0[1] = pos[3*b+1] - pos[3*a+1];
      delta0[2] = pos[3*b+2] - pos[3*a+2];
      // Apply the same cell displacement to the new relative vector
      cell_add_vec(delta0, unitcell, center);
      // Done updating delta0.
      update_delta0 = 0;
    }
    // Construct delta by adding the appropriate cell vector to delta0
    delta[0] = delta0[0];
    delta[1] = delta0[1];
    delta[2] = delta0[2];
    cell_add_vec(delta, unitcell, &(neighs[i].r0));
    // Compute the distance
    d = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
    // Store the new results;
    neighs[i].d = d;
    neighs[i].dx = delta[0];
    neighs[i].dy = delta[1];
    neighs[i].dz = delta[2];
    //neighs++;
  }
}
