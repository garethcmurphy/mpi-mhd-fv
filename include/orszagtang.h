/* $Id: initialise_uniform.h,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "out.h"

int orszagtang(TNT::Array2D<unk> mesh,
               TNT::Array2D<double> faceBx,
               TNT::Array2D<double> faceBy,
               MPI_Comm new_comm, int ndims, int *dim_size,
               double *xsize, double *ysize);

class OrszagTang {
    int nx;
    int ny;
    double gammam1i;
    double rho, bz;
    double myx;
    double myy;

    double pressure;
    double bx;
    double by;
    int coords[2];

    int myid;
    int dims[2];
    int xmax;
    int ymax;
    double vx, vy, vz, bsqr, ke;
    double centre[2];
    double a, b;
    int myaddress[2];


    unk maxt;
    unk mint;
public:
    int initial_condition(TNT::Array2D<unk> mesh,

                          TNT::Array2D<double> faceBx,
                          TNT::Array2D<double> faceBy,
                          MPI_Comm new_comm, int ndims, int *dim_size,
                          double *xsize, double *ysize
    );
};