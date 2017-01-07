/* $Id: maes.cpp,v 1.5 2006-11-16 13:48:07 gmurphy Exp $  */

/*
 * CAUTION: Needs updating
 * Has not been modified for correct addressing
 * array indices at present include boundaries
 * (need dummy variables to exclude boundaries)
 */

#include "initialise_maes.h"

int
initialise_maes(TNT::Array2D<unk> mesh,
                TNT::Array2D<double> faceBx,
                TNT::Array2D<double> faceBy,
                MPI_Comm new_comm, int ndims, int *dim_size,
                double *xsize, double *ysize) {
    double minpressure = 9e99;
    double mindensity = 9e99;
    double maxdensity = (-9e99);

    std::cout << "Initialise Magnetised Accretion-Ejection Structure " << std::
    endl;

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
    // plasma beta parameter measuring the ratio of the
    // thermal pressure to the magnetic pressure at Z = 0
    double beta = 1.0;

    double pressure = 0;
    // Disk Aspect Ratio
    double eps = 0.1;

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


    *xsize = 40;
    // Initialise the face-centred magnetic field
    for (int jj = 0; jj < ny + 1; jj++)
        for (int ii = 0; ii < nx + 1; ii++) {
            //
            int gg = std::max(ii - NGC, 0);
            float myx = 0;
            float myy = 0;
            myx = (float) myaddress[0] + gg;
            float delta_x = *xsize / ((nx - 2 * NGC) * xmax);
            //std::cout << delta_x << " " <<delta_y << std::endl;
            r = (40.0 * myx) / (dim_size[0] * nx) + delta_x * 0.5;
            h = eps * r;
            h2 = h * h;
            r2 = r * r;
            r0 = 4;
            r02 = r0 * r0;
            bz = pow(r0, 2.5) / (pow((r02 + r2), 1.25) * sqrt(beta));
            faceBx[ii][jj] = 0.0;
            faceBy[ii][jj] = bz;
            if (isnan(bz)) {
                MPI_Abort(MPI_COMM_WORLD, 33333333);
                exit(0);
            }
//      faceBy[ii][jj] = 3.0;
        }

    for (int jj = 0; jj < ny + 1; jj++) {
        faceBy[nx][jj] =
        faceBy[nx - 1][jj] =
                faceBy[nx - 2][jj];

    }

#undef PROBS_WITH_B
#ifdef PROBS_WITH_B
    for (int hh = 0; hh < ny + 1; hh++)
      for (int gg = 0; gg < nx + 1; gg++)
        {
      faceBy[hh][gg] = 3.0;
        }
#endif

    *xsize = 40;
    *ysize = 80;

    for (int jj = 2; jj < ny; jj++)
        for (int ii = 2; ii < nx; ii++) {
            // Mass Density
            int gg = ii - NGC;
            int hh = jj - NGC;
            float myx = 0;
            float myy = 0;
            myx = (float) myaddress[0] + gg;
            myy = (float) myaddress[1] + hh;
            float delta_x = *xsize / ((nx - 2 * NGC) * xmax);
            float delta_y = *ysize / ((ny - 2 * NGC) * ymax);
            //std::cout << delta_x << " " <<delta_y << std::endl;
            z = (80.0 * myy) / (dim_size[1] * ny) + delta_y * 0.5;
            r = (40.0 * myx) / (dim_size[0] * nx) + delta_x * 0.5;
            r = std::max(r, delta_x * 0.5);
            z = std::max(z, delta_y * 0.5);
            h = eps * r;
            h2 = h * h;
            z2 = z * z;
            r2 = r * r;
            /* R0 is a constant set equal to 4 in our runs.
             * This offset radius in the denominator makes the
             * density regular up to R = 0.
             */
            r0 = 4;
            r02 = r0 * r0;
            rho =
                    std::max(1e-6, (pow(r0, 1.5) / pow((r0 * r0 + r * r), 0.75))) *
                    std::pow(std::max(1e-6, (1 - 0.5 * (gammam1) * z2 / h2)), gammam1i);
            //mesh[ii][jj] _MASS= (pow(r0,1.5)/pow((r0*r0+r*r),0.75))*pow (max(1e-6,  (1- 0.5* (gammam1) * z2/h2 )), gammam1i ) ;
            rho = std::max(1e-6, rho);
            /*
            std::cout
            << rho
            << " " << ii << " " << jj
            << " " << r
            << " " <<
              std::pow (std::max (1e-6, (1 - 0.5 * (gammam1) * z2 / h2)), gammam1i)

            << std::endl;
            */


            mesh[ii][jj]_MASS = rho;
            /*
               mesh[ii][jj] _MASS =  0.2 * (
               mesh[ii][jj] _MASS  +
               mesh[ std::min(ii+1,nx-1)    ][jj] _MASS  +
               mesh[ std::max(ii-1,0)       ][jj] _MASS  +
               mesh[ii                                                      ][std::min (jj+1, ny-1) ] _MASS  +
               mesh[ii                                                      ][std ::max(jj-1,0)] _MASS
               );
             */

            // Rotation Profile
            vtheta =
                    (1 -
                     eps * eps) * sqrt(r0) * exp(-2 * z2 / h2) / (eps *
                                                                  pow((r02 + r2),
                                                                      0.25));
            mesh[ii][jj]_MOMZ = mesh[ii][jj]_MASS * vtheta;
            if (isnan(mesh[ii][jj]_MOMZ)) {
                std::cout << "z-momentum is nan" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 0);
                exit(0);
            }

            // Poloidal Velocity
            vr =
                    (-ms) * sqrt(r0) * exp(-2 * z2 / h2) / (pow((r02 + r2), 0.25));
            mesh[ii][jj]_MOMX = mesh[ii][jj]_MASS * vr;
            if (isnan(mesh[ii][jj]_MOMX)) {
                std::cout << "x-momentum is nan" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 99);
                exit(0);
            }

            vz = vr * z / r;
            mesh[ii][jj]_MOMY = mesh[ii][jj]_MOMX * z / r;
            if (isnan(mesh[ii][jj]_MOMY)) {
                std::cout << "y-momentum is nan" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 101);
                exit(0);
            }

            // Magnetic Field
            bz = pow(r0, 2.5) / (pow((r02 + r2), 1.25) * sqrt(beta));
            bz = 0.5 * (faceBy[ii][jj] + faceBy[ii][jj + 1]);
            br = 0;
            btheta = 0;

            // Pressure = rho ^ gamma
            pressure = pow(rho, PhysConsts::gamma);

            // Energy
            mesh[ii][jj]_B_X = 0.;
            mesh[ii][jj]_B_Y = bz;
            mesh[ii][jj]_B_Z = 0.;
            double vv2 = vz * vz + vr * vr + vtheta * vtheta;
            double b2 = bz * bz + br * br + btheta * btheta;
            mesh[ii][jj]_ENER =
                    (0.5 * rho * vv2) + (0.5 * b2) + pressure * gammam1i;
            if (mesh[ii][jj]_ENER < 0) {

                std::cout
                        << " Neg init energy"
                        << " (" << ii << " , " << jj << " )" << std::endl;
            }

            mesh[ii][jj].temperature = myid;

            minpressure = std::min(pressure, minpressure);
            mindensity = std::min(rho, mindensity);
            maxdensity = std::max(rho, maxdensity);
        }


    std::cout << "Min pressure " << minpressure << std::endl;
    std::cout << "Min density " << mindensity << std::endl;
    std::cout << "Max density " << maxdensity << std::endl;

    //



}
