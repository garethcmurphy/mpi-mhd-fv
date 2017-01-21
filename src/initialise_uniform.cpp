
/*
 * CAUTION: Needs updating
 * Has not been modified for correct addressing
 * array indices at present include boundaries
 * (need dummy variables to exclude boundaries)
 */

#include "initialise_uniform.h"
#include <iomanip>

int
initialise_uniform(TNT::Array2D<unk> mesh,
                   TNT::Array2D<double> faceBx,
                   TNT::Array2D<double> faceBy,
                   MPI_Comm new_comm, int ndims, int *dim_size,
                   double *xsize, double *ysize) {

    std::cout << "Initialise Orszag-Tang Problem" << std::endl;
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


    unk maxt;
    unk mint;
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

    std::cout << xmax << std::endl;
    std::cout << ymax << std::endl;
    *xsize = 1.0;
    *ysize = 1.0;
    double delta_x = *xsize / xmax;
    double delta_y = *ysize / ymax;
    // Initialise face-centered B fields
    for (int ii = 0; ii < nx + 1; ii++)
        for (int jj = 0; jj < ny + 1; jj++) {
            double myx = (double) myaddress[0] + (ii - NGC) + 0.5;
            double myy = (double) myaddress[1] + (jj - NGC) + 0.5;
            faceBx[ii][jj] = -sin(2 * PhysConsts::pi * myy / (double) ymax);
            faceBy[ii][jj] = sin(4 * PhysConsts::pi * myx / (double) xmax);
        }

    for (int jj = 2; jj < ny - 2; jj++)
        for (int ii = 2; ii < nx - 2; ii++) {
            // Mass Density

            int gg = ii - NGC;
            int hh = jj - NGC;

            bx = 0.5 * (faceBx[ii][jj] + faceBx[ii + 1][jj]);
            by = 0.5 * (faceBy[ii][jj] + faceBy[ii][jj + 1]);
            bz = 0.0;
            rho = 5 / (12 * PhysConsts::pi);
            rho = 25. / 9.;
            pressure = 25 / (36 * PhysConsts::pi);
            pressure = 5. / 3.;
            //rho=0.2210485;
            //pressure = 0.13269;
            double myx = (float) myaddress[0] + (ii - NGC) + 0.5;
            double myy = (float) myaddress[1] + (jj - NGC) + 0.5;


//	std::cout << myx << " "<< myy << std::endl;

            double vx = -sin(2.0 * PhysConsts::pi * myy / ymax);
            //std::cout << std::setiosflags(std::ios::scientific);
            //std::cout << rho *vx << " "<< myy << " "  << ymax << std::endl;
            double vy = sin(2.0 * PhysConsts::pi * myx / xmax);
            double vz = 0.0;

            double bsqr = 0.5 * (bx * bx + by * by + bz * bz);
            double ke = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
            mesh[ii][jj]_MASS = rho;
            mesh[ii][jj]_MOMX = rho * vx;
            mesh[ii][jj]_MOMY = rho * vy;
            mesh[ii][jj]_MOMZ = 0.0;
            mesh[ii][jj]_ENER = pressure * gammam1i + bsqr + ke;
            mesh[ii][jj]_B_X = bx;
            mesh[ii][jj]_B_Y = by;
            mesh[ii][jj]_B_Z = bz;
            mesh[ii][jj].temperature = myid;


            float centre[2];
            centre[0] = xmax * (nx - 2 * NGC) / 2 + 40;
            centre[1] = ymax * (ny - 2 * NGC) / 2 + 40;

        }

    for (int jj = 0; jj < ny; jj++) {
        mesh[0][jj] = mesh[nx - 4][jj];
        mesh[1][jj] = mesh[nx - 3][jj];
        mesh[nx - 1][jj] = mesh[3][jj];
        mesh[nx - 2][jj] = mesh[2][jj];
    }
    for (int ii = 0; ii < nx; ii++) {
        mesh[ii][0] = mesh[ii][4];
        mesh[ii][1] = mesh[ii][3];
        mesh[ii][ny - 1] = mesh[ii][3];
        mesh[ii][ny - 1] = mesh[ii][2];
    }
    for (int jj = 2; jj < ny - 2; jj++)
        for (int ii = 2; ii < nx - 2; ii++) {

            double myx = (float) myaddress[0] + (ii - NGC) + 0.5;
            double myy = (float) myaddress[1] + (jj - NGC) + 0.5;
            double a = mesh[ii][jj]_MOMX;
            double b = -mesh[nx - ii - 1][ny - jj - 1]_MOMX;

            if (a - b > 1e-10) {

                std::cout
                        << ii << " " << jj
                        //<< " " << nx-ii-1 << " " << ny-jj-1
                        //<< " " << myx << " " << myy
                        //<< " " << a << " " << b
                        << " " << a - b
                        << std::endl;
            }

            for (int qq = 0; qq < 9; ++qq) {
                maxt.array[qq] = maxt.array[qq] > mesh[ii][jj].array[qq] ? maxt.array[qq] : mesh[ii][jj].array[qq];
                mint.array[qq] = mint.array[qq] < mesh[ii][jj].array[qq] ? mint.array[qq] : mesh[ii][jj].array[qq];
            }


        }


    for (int qq = 0; qq < 9; ++qq) {
        std::cout <<
                  " Max, min = " <<
                  maxt.array[qq]
                  << " , " <<
                  mint.array[qq]
                  << std::endl;
    }
    return 0;

}
