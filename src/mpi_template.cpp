#include <mpi.h>
int
main (int argc, char *argv[])
{
  int numprocs = 0;
  int myid = 0;
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);

  MPI_Finalize ();

  return 0;
}
