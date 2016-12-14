/* $Id: initialise_jet.cpp,v 1.5 2006-11-16 13:48:06 gmurphy Exp $  */

/*
 */

#include "initialise_jet.h"

int
initialise_jet(TNT::Array2D<unk> mesh,
               TNT::Array2D<double> faceBx,
               TNT::Array2D<double> faceBy,
               MPI_Comm new_comm, int ndims, int *dim_size,
               double *xsize, double *ysize) {

    std::cout << "Initialise Jet Problem" << std::endl;
    const int nx = mesh.dim1();
    const int ny = mesh.dim2();
    const double gammam1i = 1 / (PhysConsts::gamma - 1);
    const double gammam1 = (PhysConsts::gamma - 1);
    double r0, r, z, h;
    double r2, z2, h2;
    double r02;
    // ms is a parameter smaller than unity
    // which then ensures an initial subsonic poloidal
    // inflow
    double ms = 0.3;
    double rho, vtheta, vz, vr, bz, br, btheta;
    // plasma beta parameter measuring the ratio of the thermal pressure to the magnetic pressure at Z = 0
    double beta = 1.0;

    double pressure = 0;
    // Disk Aspect Ratio
    double eps = 0.1;
    double bx = 10;
    double by = 10;


    int coords[] = {0, 0};

    int myid = 99;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    MPI_Cart_coords(new_comm, myid, ndims, coords);
    int dims[] = {0, 0};
    MPI_Cartdim_get(new_comm, dims);

    int xmax = 0;
    int ymax = 0;
    MPI_Allreduce(coords, &xmax, 1, MPI_INT, MPI_MAX, new_comm);
    MPI_Allreduce(coords + 1, &ymax, 1, MPI_INT, MPI_MAX, new_comm);
    xmax++;
    ymax++;

    std::cout << "Proc" << myid << " " << xmax << " " << ymax << std::endl;
    int myaddress[] = {0, 0};
    myaddress[0] = (nx - 2 * NGC) * coords[0];
    myaddress[1] = (ny - 2 * NGC) * coords[1];

    xmax *= (nx - 2 * NGC);
    ymax *= (ny - 2 * NGC);

    *xsize = 100.0;
    *ysize = 200.0;
    // Initialise face-centered B fields
    for (int ii = 0; ii < nx + 1; ii++)
        for (int jj = 0; jj < ny + 1; jj++) {
            int gg = ii - NGC;
            int hh = jj - NGC;
            double myx = (float) myaddress[0] + (ii - NGC);
            double myy = (float) myaddress[1] + (jj - NGC);
            //faceBx[ii][jj] =  -sin( 2*PhysConsts::pi * myy / ymax) ;
            //faceBy[ii][jj] =  sin( 4*PhysConsts::pi * myx / xmax);
            faceBx[ii][jj] = 0.0;
            faceBy[ii][jj] = 1.0;
        }

    for (int jj = 0; jj < ny; jj++)
        for (int ii = 0; ii < nx; ii++) {
            // Mass Density
            bx = 0.5 * (faceBx[ii][jj] + faceBx[ii + 1][jj]);
            by = 0.5 * (faceBy[ii][jj] + faceBy[ii][jj + 1]);
            bz = 0.0;
            rho = 0.25;
            double myx = (float) myaddress[0] + (ii - NGC);
            double myy = (float) myaddress[1] + (jj - NGC);
            double vx = 0.0;
            double vy = 0.0;
            double vz = 0.0;
            if (0 && myx < 40 && myy < 40) {
                vy = 10.0;
                if (myx > 30) {
                    vy = (40 - 30) * (myx - 40) / (0 - 10);
                }
            }
            //if (  myy<40) vy=10.0;

            double bsqr = 0.5 * (bx * bx + by * by + bz * bz);
            double ke = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
            double rmax2 = nx * nx;
            pressure = 0.5 + (1. - myx * myx / rmax2);
            pressure = 0.0004;
            mesh[ii][jj]_MASS = rho;
            mesh[ii][jj]_MOMX = rho * vx;
            mesh[ii][jj]_MOMY = rho * vy;
            mesh[ii][jj]_MOMZ = rho * vz;
            mesh[ii][jj]_ENER = pressure * gammam1i + bsqr + ke;
            mesh[ii][jj]_B_X = bx;
            mesh[ii][jj]_B_Y = by;
            mesh[ii][jj]_B_Z = bz * myx / nx;
            mesh[ii][jj].temperature = myid;
        }


}
