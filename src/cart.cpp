/* $Id: cart.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include <mpi.h>
#include <iostream>
using namespace std;
int
main (int argc, char *argv[])
{
  int numprocs = 0;
  int myid = 0;

  MPI_Comm old_comm, new_comm;
  int ndims, reorder, ierr;
  int dim_size[2], periods[2];
  int ii = 0, jj = 0;
  int coords[] = { 0, 0 };
  int rank = 0;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);

  old_comm = MPI_COMM_WORLD;
  ndims = 2;
  dim_size[0] = 3;
  dim_size[1] = 2;
  periods[0] = 1;
  periods[1] = 0;
  reorder = 1;

  MPI_Cart_create (old_comm, ndims, dim_size, periods, reorder, &new_comm);
  for (ii = 0; ii < 3; ii++)
    {
      for (jj = 0; jj < 3; jj++)
	{
	  coords[0] = ii;
	  coords[1] = jj;
	  MPI_Cart_rank (new_comm, coords, &rank);
	  cout << " " << coords[1] << " " << rank << endl;
	}
    }

  MPI_Finalize ();

  return 0;
}
