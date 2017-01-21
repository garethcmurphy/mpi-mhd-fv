#include "out.h"
int initialise_blast (TNT::Array2D < unk > mesh,
		      TNT::Array2D < double >faceBx,
		      TNT::Array2D < double >faceBy,
		      MPI_Comm new_comm, int ndims, int *dim_size,
		      double *xsize, double *ysize);


class Blast {
    int nx;
    int ny;
    double gammam1i;
    double pressure;
    double bx;
    double by;
    double sqr4piei;
    double myx;
    double myy;
    double rh;
    double vx;
    double vy;
    double vz;
    double centre[2];
    int coords[2];
    int dims[2];
    int myid;
    int xmax;
    int ymax;

public:
    int initial_condition(TNT::Array2D<unk> mesh,
                          TNT::Array2D<double> faceBx,
                          TNT::Array2D<double> faceBy,
                          MPI_Comm new_comm, int ndims);
};

