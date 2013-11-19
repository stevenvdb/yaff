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
#include "mpi.h"
#include "fcs.h"

double compute_ewald_reci_scafacos(double *pos, long natom, double *charges,
                          cell_type* cell, double alpha, long *gmax,
                          double gcut, double *gpos, double *work,
                          double* vtens) {
    int comm_size;
    int comm_rank;
  fcs_int periodicity[3];
  fcs_float *box_a;
  fcs_float *box_b;
  fcs_float *box_c;
  fcs_float *offset;
  fcs_float *particles;
  fcs_float *sfccharges;
  fcs_float *velocities;
  fcs_float *masses;
  fcs_float *field;
  fcs_float *potentials;
  fcs_int local_particles;


    int i;
    MPI_Comm communicator = MPI_COMM_WORLD;
    MPI_Init (NULL, NULL);
    MPI_Comm_size (communicator, &comm_size);
    printf("\nEWALD COMM SIZE: %d\n", comm_size);

  // Initialization
  char *method;
  //method = "ewald";
  method = "ewald";
  FCS handle = NULL;
  FCSResult result = NULL;
  result = fcs_init (&handle, method, communicator);
  MPI_Comm_rank (communicator, &comm_rank);
  if (comm_rank == 0){
    printf ("Initializing FCS...\n\n");
    fcs_result_print_result (result);
    printf("\n");
  }
  fcs_result_destroy (result);
  MPI_Barrier(communicator);

  // Setup
  box_a = (fcs_float *) malloc (3 * sizeof (fcs_float));
  box_b = (fcs_float *) malloc (3 * sizeof (fcs_float));
  box_c = (fcs_float *) malloc (3 * sizeof (fcs_float));
  offset = (fcs_float *) malloc (3 * sizeof (fcs_float));
  for (i=0;i<3;i++){
     box_a[i] = (*cell).rvecs[i];
     box_b[i] = (*cell).rvecs[3+i];
     box_c[i] = (*cell).rvecs[6+i];
     offset[i] = 0.0;
     periodicity[i] = 1;
  }
  fcs_set_common (handle, 1, box_a, box_b, box_c, offset, periodicity, natom);

  /* dividing particles up */
  if (comm_rank == comm_size - 1)
    local_particles =
      natom - (comm_size -
                 1) * (natom / comm_size);
  else
    local_particles = natom / comm_size;



  // Tuning
  particles = (fcs_float *) malloc (3 * sizeof (fcs_float) * local_particles);
  sfccharges = (fcs_float *) malloc (3 * sizeof (fcs_float) * local_particles);
  velocities =
    (fcs_float *) malloc (3 * sizeof (fcs_float) * local_particles);
  masses = (fcs_float *) malloc (sizeof (fcs_float) * local_particles);
  field = (fcs_float *) malloc (3 * sizeof (fcs_float) * local_particles);
  potentials = (fcs_float *) malloc (sizeof (fcs_float) * local_particles);

  int j;
  j = 0;
  for (i = 0; i < natom; ++i)
    {
      if (((i >= comm_rank * local_particles)
       || i >= natom - local_particles)
      && i < (comm_rank + 1) * local_particles)
    {
      particles[3*j+0] = pos[3*i+0];
      particles[3*j+1] = pos[3*i+1];
      particles[3*j+2] = pos[3*i+2];
      sfccharges[j] = charges[i];
      velocities[3*j+0] = 0.0;
      velocities[3*j+1] = 0.0;
      velocities[3*j+2] = 0.0;
      masses[j] = 1.0;
      j++;
    }
    }

  for (i = 0; i < comm_size; ++i)
     {
     if (comm_rank == i)
     {
     for (j = 0; j < local_particles; ++j)
     printf("%d %e %e %e %e\n", comm_rank, *(particles+j*3), *(particles+j*3+1), *(particles+j*3+2), *(charges+j));
     }
     MPI_Barrier(communicator);
     }

  result = fcs_ewald_set_alpha( handle, alpha);
  result = fcs_ewald_set_kmax( handle, gmax[0]);
  //result = fcs_ewald_set_r_cut( handle, 10.0);
  //result = fcs_p3m_set_alpha( handle, alpha);
  if (comm_rank == 0) {
     fprintf (stderr, "Setting parameters for FCS...\n\n");
     fcs_result_print_result (result);
     fprintf(stderr, "\n");
  }
  fcs_result_destroy (result);

  if (comm_rank == 0){
    printf ("\n");
    printf ("Parameters... \n");
    fcs_print_parameters (handle);
  }
  MPI_Barrier(communicator);

  // Tuning step
  MPI_Barrier(communicator);

  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Tuning... \n");
  result =
    fcs_tune (handle, local_particles, particles, sfccharges);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);

  MPI_Barrier(communicator);


  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Parameters... \n");
  if (comm_rank == 0){
    fcs_print_parameters (handle);
    printf ("\n");
  }

   MPI_Barrier(communicator);

   // Running Step
  if (comm_rank == 0){
    printf ("\n");
    printf ("Running... \n");}


   result = fcs_run (handle, local_particles, particles,
            sfccharges, field, potentials);

  if (comm_rank == 0){
    printf ("\n");
    printf ("Done Running... \n");}



  double E = 0.0;
  for (i = 0; i < comm_size; ++i){
     if (i == comm_rank)
     {
          for (j = 0; j < local_particles; ++j){
          printf("%d %e %e %e %e\n", comm_rank, *(field+j*3), *(field+j*3+1), *(field+j*3+2), *(potentials+j));
          E += (*(potentials+j)) * (*(charges+j));
          gpos[3*j+0] -= (*(field+j*3+0)) * (*(charges+j));
          gpos[3*j+1] -= (*(field+j*3+1)) * (*(charges+j));
          gpos[3*j+2] -= (*(field+j*3+2)) * (*(charges+j));
          }
     }
     MPI_Barrier(communicator);
     }



   // Destroy step
  if (comm_rank == 0){
    printf ("\n");
    printf ("Destroying... \n");}
  fcs_result_destroy (result);



  /*
  free (box_a);
  free (box_b);
  free (box_c);
  free (offset);

  fcs_destroy (handle);

  free (particles);
  free (fsccharges);
  free (velocities);
  free (field);
  free (potentials);
  free (masses);

  //MPI_Finalize ();

  */


    return 0.5*E;
}




def compute_p3m_scafacos(np.ndarray[double, ndim=2] pos,
                       np.ndarray[double, ndim=1] charges,
                       Cell unitcell,
                       np.ndarray[double, ndim=2] gpos,
                       np.ndarray[double, ndim=1] work,
                       np.ndarray[double, ndim=2] vtens):
    '''Compute the reciprocal interaction term in the Ewald summation scheme

       **Arguments:**

       pos
            The atomic positions. numpy array with shape (natom,3).

       charges
            The atomic charges. numpy array with shape (natom,).

       unitcell
            An instance of the ``Cell`` class that describes the periodic
            boundary conditions.

       gpos
            If not set to None, the Cartesian gradient of the energy is
            stored in this array. numpy array with shape (natom, 3).

       work
            If gpos is given, this work array must also be present. Its
            contents will be overwritten. numpy array with shape (2*natom,).

       vtens
            If not set to None, the virial tensor is computed and stored in
            this array. numpy array with shape (3, 3).
    '''
    cdef double *my_gpos
    cdef double *my_work
    cdef double *my_vtens

    assert pos.flags['C_CONTIGUOUS']
    assert pos.shape[1] == 3
    assert charges.flags['C_CONTIGUOUS']
    assert charges.shape[0] == pos.shape[0]
    assert unitcell.nvec == 3

    if gpos is None:
        my_gpos = NULL
        my_work = NULL
    else:
        assert gpos.flags['C_CONTIGUOUS']
        assert gpos.shape[1] == 3
        assert gpos.shape[0] == pos.shape[0]
        assert work.flags['C_CONTIGUOUS']
        assert gpos.shape[0]*2 == work.shape[0]
        my_gpos = <double*>gpos.data
        my_work = <double*>work.data

    if vtens is None:
        my_vtens = NULL
    else:
        assert vtens.flags['C_CONTIGUOUS']
        assert vtens.shape[0] == 3
        assert vtens.shape[1] == 3
        my_vtens = <double*>vtens.data
    p3m = scafacos.scafacos("p3m")
    print "HELLO"
    p3m.box = unitcell.rvecs
    print "Setting periodicity"
    p3m.periodicity = (True, True, True)
    print "Setting nead field flag"
    p3m.near_field_flag = True
    print "Setting number of parricls"
    positions = np.array([[0.5, 1.0, 1.0], [1.5, 1.0, 1.0],[0.2,0.2,0.2],[0.9,0.9,0.9]], dtype='double')
    charges = np.array([1.0, -1.0,0.5,-0.5], dtype='double')
    p3m.total_particles = charges.shape[0]
    print "Tuning"
    p3m.tune(pos, charges)
    print "Computing"
    fields, potentials = p3m(pos, charges)
    print fields
    print potentials
    return 0.0
