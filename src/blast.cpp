
/*
 * CAUTION: Needs updating
 * Has not been modified for correct addressing
 * array indices at present include boundaries
 * (need dummy variables to exclude boundaries)
 */

#include "initialise_blast.h"

int
initialise_blast(TNT::Array2D<unk> mesh,
                 TNT::Array2D<double> faceBx,
                 TNT::Array2D<double> faceBy,
                 MPI_Comm new_comm, int ndims, int *dim_size,
                 double *xsize, double *ysize) {

    *xsize = 1.0;
    *ysize = 1.0;
    std::cout << "Initialise Blast Problem" << std::endl;
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
    double sqr4piei = 0.5 / std::sqrt(PhysConsts::pi);
    bx = bx * sqr4piei;
    by = by * sqr4piei;


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


    // Initialise face-centered B fields
    for (int ii = 0; ii < nx + 1; ii++)
        for (int jj = 0; jj < ny + 1; jj++) {
            faceBx[ii][jj] = bx;
            faceBy[ii][jj] = by;
        }

    for (int jj = 0; jj < ny; jj++)
        for (int ii = 0; ii < nx; ii++) {
            // Mass Density

            int gg = ii - NGC;
            int hh = jj - NGC;

            bx = 0.5 * (faceBx[ii][jj] + faceBx[ii + 1][jj]);
            by = 0.5 * (faceBy[ii][jj] + faceBy[ii][jj + 1]);
            pressure = 0.1;

            mesh[ii][jj]_MASS = 1.0;
            mesh[ii][jj]_MOMX = 0.0;
            mesh[ii][jj]_MOMY = 0.0;
            mesh[ii][jj]_MOMZ = 0.0;
            mesh[ii][jj]_ENER = pressure * gammam1i + 0.5 * (bx * bx + by * by);
            mesh[ii][jj]_B_X = bx;
            mesh[ii][jj]_B_Y = by;
            mesh[ii][jj]_B_Z = 0.0;
            mesh[ii][jj].temperature = myid;

            float myx = 0;
            float myy = 0;

            float centre[2];
            centre[0] = xmax * (nx - 2 * NGC) / 2 + 10;
            centre[1] = ymax * (ny - 2 * NGC) / 2 + 10;


            myx = (float) myaddress[0] + gg - centre[0];
            myy = (float) myaddress[1] + hh - centre[1];

            double dist = 0;
            dist = sqrt((myx * myx) + (myy * myy));
/* myx = (float) myaddress[0] + gg ;
	myy = (float) myaddress[1] + hh ;
	dist= sqrt ((myx * myx));
	*/

            if (dist < 0.1 * nx)
//          if (sqrt ( (myy * myy)) < 50.0)
            {
                double rh = 1.0;
                double vx = 0;
                double vy = 0.0;
                double vz = 0.0;

                mesh[ii][jj]_MASS = rh;
                mesh[ii][jj]_MOMX = vx;
                mesh[ii][jj]_MOMY = vy;
                mesh[ii][jj]_MOMZ = vz;
                mesh[ii][jj]_B_X = bx;
                mesh[ii][jj]_B_Y = by;
                mesh[ii][jj]_B_Z = 0.0;
                mesh[ii][jj]_ENER = 10.0 * gammam1i
                                    + 0.5 * rh * (vx * vx + vy * vy) + 0.5 * (bx * bx + by * by);

            }


        }


}
