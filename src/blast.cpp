
#include "initialise_blast.h"


int Blast::initial_condition(TNT::Array2D<unk> mesh,
                             TNT::Array2D<double> faceBx,
                             TNT::Array2D<double> faceBy,
                             MPI_Comm new_comm, int ndims) {

    nx = mesh.dim1();
    ny = mesh.dim2();
    gammam1i = 1 / (PhysConsts::gamma - 1);
    pressure = 0;
    bx = 10;
    by = 10;
    sqr4piei = 0.5 / std::sqrt(PhysConsts::pi);
    myx = 0;
    myy = 0;
    rh = 1.0;
    vx = 0;
    vy = 0.0;
    vz = 0.0;
    myid = 99;
    xmax = 0;
    ymax = 0;

    bx = bx * sqr4piei;
    by = by * sqr4piei;


    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    MPI_Cart_coords(new_comm, myid, ndims, coords);
    MPI_Cartdim_get(new_comm, dims);

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

            centre[0] = xmax * (nx - 2 * NGC) / 2 + 10;
            centre[1] = ymax * (ny - 2 * NGC) / 2 + 10;


            myx = (float) myaddress[0] + gg - centre[0];
            myy = (float) myaddress[1] + hh - centre[1];

            double dist = 0;
            dist = sqrt((myx * myx) + (myy * myy));

            if (dist < 0.1 * nx) {

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


    return 0;
};

