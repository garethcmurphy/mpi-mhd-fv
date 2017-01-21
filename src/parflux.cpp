
#include "mpi.h"
#include "tnt.h"
#include "out.h"
#include "parflux.h"
#include "riemann.h"
#include "hlld.h"
#include "problem.h"
#include "physics.h"
#include "ctop.h"

int
errcheck(int rc, int *myaddress, int ii, int jj) {
    if (rc == 1) {
        std::
        cout << __FUNCTION__ << ": location on proc " << ii << "," << jj <<
             std::endl;
        double myx = (double) myaddress[0] + ii;
        double myy = (double) myaddress[1] + jj;
        std::
        cout << __FUNCTION__ << ": global location " << myx << "," << myy <<
             std::endl;
        MPI_Abort(MPI_COMM_WORLD, 50505050);
        exit(0);
    }
    return 0;
}

int
parflux(TNT::Array2D<unk> mesh,
        TNT::Array2D<flux> fx,
        TNT::Array2D<flux> fy,
        TNT::Array2D<double> xjx,
        TNT::Array2D<double> xjy,
        TNT::Array2D<double> xjz,
        TNT::Array2D<double> yjx,
        TNT::Array2D<double> yjy,
        TNT::Array2D<double> yjz,
        TNT::Array2D<double> faceEx,
        TNT::Array2D<double> faceEy,
        TNT::Array2D<double> faceBx,
        TNT::Array2D<double> faceBy,
        TNT::Array1D<double> SoundSpeedMidplane,
        double del,
        double delta_x, double delta_y, int *myaddress,
        int tstep,
        MPI_Comm Cart_comm) {

    double eps = 0.1;
    double gammam1 = PhysConsts::gamma - 1.;
    double gammam1i = 1. / gammam1;

    int nx = mesh.dim1();
    int ny = mesh.dim2();
    // Current arrays

    int myid = 99999;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    int ndims = 2;
    int coords[2];


    int remain_dims[] = {0, 1};
    MPI_Comm Sub_comm;
    MPI_Cart_sub(Cart_comm, remain_dims, &Sub_comm);
    int lowest = 99999;
    coords[1] = 1;
    MPI_Cart_coords(Cart_comm, myid, ndims, coords);
    MPI_Cart_rank(Cart_comm, coords, &lowest);
    // Check my position in the Cartesian grid
    // Then bcast my id to the Sub_comm
    // Then all print out the id
    // Grab the coords of the processor
    double tmpbuf;
    tmpbuf = myid;
    TNT::Array1D<double> pressure_buf(nx);
    MPI_Bcast(&tmpbuf, 1, MPI_DOUBLE, 0, Sub_comm);

    SoundSpeedMidplane = 0;
    for (int i = 20; i < nx; ++i) {
        double rr = mesh[i][2]_MASS;
        double px = mesh[i][2]_MOMX;
        double py = mesh[i][2]_MOMY;
        double pz = mesh[i][2]_MOMZ;
        double et = mesh[i][2]_ENER;
        double bx = mesh[i][2]_B_X;
        double by = mesh[i][2]_B_Y;
        double bz = mesh[i][2]_B_Z;

        double ri = 1.0 / rr;
        double vx = px * ri;
        double vy = py * ri;
        double vz = pz * ri;
        double ke = 0.5 * rr * (vx * vx + vy * vy + vz * vz);
        double b2 = 0.5 * (bx * bx + by * by + bz * bz);
        double p = et - ke - b2;
        double pmin = 1e-20;
        p = std::max(p, pmin);
        p = p * gammam1;
        double cs = std::sqrt(PhysConsts::gamma * p * ri);
        SoundSpeedMidplane[i] = cs;
    }
    int rc = MPI_Bcast(&SoundSpeedMidplane[0], nx, MPI_DOUBLE, 0, Sub_comm);
    if (rc != 0) {
        std::cerr << "Bcast fail" << std::endl;
    }

#ifdef DEBUG_SUB_COMM
    std::cout
      << "myid = " << myid
      << "coords = (" << coords[0] << "," << coords[1] << ")"
      << " lowest = " << tmpbuf << std::endl;
#endif


    int loc[2];
    double rz[2];
    double leftprim[8];
    double rightprim[8];

#ifdef RESISTIVE
    for (int ii = 1; ii < nx - 1; ii++)
      {
        for (int jj = 1; jj < ny - 1; jj++)
      {

        double rr = mesh[ii][jj] _MASS;
        if (rr < 0)
          {
            std::cout
          << __FUNCTION__ << ": "
          << "Negative density at"
          << " (" << ii << " , " << jj << " )" << std::endl;
            MPI_Abort (MPI_COMM_WORLD, 1234);
            exit (0);
          }
        double px = mesh[ii][jj] _MOMX;
        double py = mesh[ii][jj] _MOMY;
        double pz = mesh[ii][jj] _MOMZ;
        double et = mesh[ii][jj] _ENER;
        double bx = mesh[ii][jj] _B_X;
        double by = mesh[ii][jj] _B_Y;
        double bz = mesh[ii][jj] _B_Z;

        double ri = 1.0 / rr;
        double vx = px * ri;
        double vy = py * ri;
        double vz = pz * ri;
        double ke = 0.5 * rr * (vx * vx + vy * vy + vz * vz);
        double b2 = 0.5 * (bx * bx + by * by + bz * bz);
        double vdotb = (vx * bx + vy * by + vz * bz);
        double p = et - ke - b2;
        double ptot = et - ke;
        p = p * gammam1;


        if (isnan (p))
          {
            std::cout
          << __FUNCTION__
          << " Pressure is nan "
          << " ("
          << ii << "," << jj
          << ") bx " << bx
          << " by " << by
          << " bz " << bz
          << " b2 " << b2
          << " et " << et
          << std::endl;
            MPI_Abort (MPI_COMM_WORLD, 44444);
            exit (0);
          }

        // Work out current at each cell face and add to energy flux
        // xflux current
        double idx = 1.0 / (delta_x);
        double idy = 1.0 / (delta_y);


        xjx[ii][jj] =
          0.5 * idy * ((mesh[ii][jj + 1] _B_Z + mesh[ii - 1][jj + 1] _B_Z) -
               (mesh[ii][jj - 1] _B_Z + mesh[ii - 1][jj - 1] _B_Z));

        if (isnan (xjx[ii][jj]))
          {
            std::cout
          << "xjx"
          << " " << delta_y
          << " " << idy
          << " " << delta_x
          << " " << idx
          << " " << mesh[ii][jj + 1] _B_Z
          << " " << mesh[ii - 1][jj + 1] _B_Z << std::endl;
            MPI_Abort (MPI_COMM_WORLD, 666666);
            exit (0);
          }

        xjy[ii][jj] =
          (-idx) * ((mesh[ii][jj] _B_Z - mesh[ii - 1][jj] _B_Z));
        xjz[ii][jj] =
          (idx) * ((mesh[ii][jj] _B_Y - mesh[ii - 1][jj] _B_Y)) -
          0.5 * idy * ((mesh[ii][jj + 1] _B_X + mesh[ii - 1][jj + 1] _B_X) -
               (mesh[ii][jj - 1] _B_X + mesh[ii - 1][jj - 1] _B_X));

        // yflux current
        yjx[ii][jj] = idy * (mesh[ii][jj] _B_Z - mesh[ii][jj - 1] _B_Z);
        yjy[ii][jj] =
          0.5 * (-idx) *
          ((mesh[ii + 1][jj] _B_Z + mesh[ii + 1][jj - 1] _B_Z) -
           (mesh[ii - 1][jj] _B_Z + mesh[ii - 1][jj - 1] _B_Z));
        yjz[ii][jj] =
          0.5 * (idx) *
          ((mesh[ii + 1][jj] _B_Y + mesh[ii + 1][jj - 1] _B_Y) -
           (mesh[ii - 1][jj] _B_Y + mesh[ii - 1][jj - 1] _B_Y)) -
          idy * (mesh[ii][jj] _B_X - mesh[ii][jj - 1] _B_X);

      }
      }
#endif

    // Sweep for x fluxes
    for (int ii = 2; ii < nx - 1; ii++) {
        for (int jj = 2; jj < ny - 2; jj++) {
            double lstate[8];
            double rstate[8];
            double iflux[8];
            double Res_state[8];
            int timestep = tstep;
            double unused = 0;


            int rc = 0;

#ifdef ROE
            rc =
              riemann (mesh[ii - 1][jj].array, mesh[ii][jj].array,
                   fx[ii][jj].array, Res_state, timestep, &unused, 1);
            errcheck (rc, myaddress, ii, jj);
#endif

#ifdef HLLD

            loc[0] = (double) myaddress[0] + ii - 1;
            loc[1] = (double) myaddress[1] + jj;

            rz[0] = loc[0] * delta_x;
            rz[1] = loc[1] * delta_y;

            rc = ctop(mesh[ii - 1][jj].array, leftprim, rz);
            errcheck(rc, myaddress, ii - 1, jj);

            loc[0] = (double) myaddress[0] + ii;
            loc[1] = (double) myaddress[1] + jj;
            rz[0] = loc[0] * delta_x;
            rz[1] = loc[1] * delta_y;

            rc = ctop(mesh[ii][jj].array, rightprim, rz);
            errcheck(rc, myaddress, ii, jj);

            //leftprim[4]= pow( leftprim[0] , PhysConsts::gamma);
            //rightprim[4]= pow( rightprim[0] , PhysConsts::gamma);
#ifdef STAGGER_MESH
            leftprim[5] = faceBx[ii][jj];
            rightprim[5] = faceBx[ii][jj];
#endif
            rc = hlld(leftprim, rightprim, fx[ii][jj].array, Res_state, timestep, &unused, 1, loc);
#ifdef LAXFRIEDRICHS
            if (fabs (leftprim[5]) > 1e10
                ||  fabs (leftprim[6]) > 1e10
                ||  fabs (rightprim[6]) > 1e10)
            {
                std::cout
                << "At " << ii << "," << jj
                << " faceBx = " << leftprim[5]
                << " By = " << leftprim[6]
                << " By = " << rightprim[6]
                << std::endl;
            }

              rc = lf (leftprim, rightprim, fx[ii][jj].array, Res_state, timestep, &unused, 1, loc);
#endif
            errcheck(rc, myaddress, ii, jj);
#endif
        }
    }
    std::cout
            << " 5"
            << " 5"
            << " " << mesh[4][5]_MASS
            << " " << mesh[5][5]_MASS
            << " " << faceBx[5][5]
            << std::endl;

    // Sweep for y fluxes
    for (int ii = 2; ii < nx - 2; ii++) {
        for (int jj = 2; jj < ny - 1; jj++) {

            double lstate[8];
            double rstate[8];
            double iflux[8];
            double Res_state[8];
            int timestep = tstep;
            double unused = 0;
            lstate[0] = mesh[ii][jj - 1].array[0];
            lstate[1] = mesh[ii][jj - 1].array[2];
            lstate[2] = mesh[ii][jj - 1].array[3];
            lstate[3] = mesh[ii][jj - 1].array[1];
            lstate[4] = mesh[ii][jj - 1].array[4];
            lstate[5] = mesh[ii][jj - 1].array[6];
            lstate[6] = mesh[ii][jj - 1].array[7];
            lstate[7] = mesh[ii][jj - 1].array[5];

            rstate[0] = mesh[ii][jj].array[0];
            rstate[1] = mesh[ii][jj].array[2];
            rstate[2] = mesh[ii][jj].array[3];
            rstate[3] = mesh[ii][jj].array[1];
            rstate[4] = mesh[ii][jj].array[4];
            rstate[5] = mesh[ii][jj].array[6];
            rstate[6] = mesh[ii][jj].array[7];
            rstate[7] = mesh[ii][jj].array[5];

#ifdef HLLD
            loc[0] = (double) myaddress[0] + ii;
            loc[1] = (double) myaddress[1] + jj - 1;
            rz[0] = loc[0] * delta_x;
            rz[1] = loc[1] * delta_y;
            rc = ctop(lstate, leftprim, rz);
            errcheck(rc, myaddress, ii, jj - 1);

            loc[0] = (double) myaddress[0] + ii;
            loc[1] = (double) myaddress[1] + jj;
            rz[0] = loc[0] * delta_x;
            rz[1] = loc[1] * delta_y;
            rc = ctop(rstate, rightprim, rz);
            errcheck(rc, myaddress, ii, jj);
            //leftprim[4]= pow( leftprim[0] , PhysConsts::gamma);
            //rightprim[4]= pow( rightprim[0] , PhysConsts::gamma);
#ifdef STAGGER_MESH
            leftprim[5] = faceBy[ii][jj];
            rightprim[5] = faceBy[ii][jj];
#endif
            rc = hlld(leftprim, rightprim, iflux, Res_state, timestep, &unused, 2, loc);
#ifdef LAXFRIEDRICHS
            if (fabs (leftprim[5]) > 1e10)
            {
                std::cout
                << "At " << ii << "," << jj
                << " faceBx = " << leftprim[5]
                << std::endl;
            }
              rc = lf (leftprim, rightprim, iflux, Res_state, timestep, &unused, 1, loc);
#endif
#endif
#ifdef ROE
            rc =
              riemann (lstate, rstate, iflux, Res_state, timestep, &unused, 1);
            errcheck (rc, myaddress, ii, jj);
#endif
            fy[ii][jj].array[0] = iflux[0];
            fy[ii][jj].array[2] = iflux[1];
            fy[ii][jj].array[3] = iflux[2];
            fy[ii][jj].array[1] = iflux[3];
            fy[ii][jj].array[4] = iflux[4];
            fy[ii][jj].array[6] = iflux[5];
            fy[ii][jj].array[7] = iflux[6];
            fy[ii][jj].array[5] = iflux[7];
        }
    }


#undef LAPIDUS_ARTIFICIAL_VISCOSITY
#ifdef LAPIDUS_ARTIFICIAL_VISCOSITY

    for (int ii = 1; ii < nx - 1; ii++)
      {
        for (int jj = 1; jj < ny - 1; jj++)
      {

        double u1 = mesh[ii + 1][jj] _MOMX / mesh[ii + 1][jj] _MASS;
        double u2 = mesh[ii - 1][jj] _MOMX / mesh[ii - 1][jj] _MASS;
        double v1 = mesh[ii][jj + 1] _MOMY / mesh[ii][jj + 1] _MASS;
        double v2 = mesh[ii][jj - 1] _MOMY / mesh[ii][jj - 1] _MASS;
        double divv = 0.5 * (u1 - u2 + v1 - v2);
        double delu = 0;


        for (int hh = 0; hh < 8; hh++)
          {
            delu = (mesh[ii][jj].array[hh] - mesh[ii - 1][jj].array[hh]);
            fx[ii][jj].array[hh] += del * 0.1 * min (0.0, divv) * delu;
            delu = (mesh[ii][jj].array[hh] - mesh[ii][jj - 1].array[hh]);
            fy[ii][jj].array[hh] += del * 0.1 * min (0.0, divv) * delu;
          }

      }
      }
#endif


#ifdef RESISTIVE
    for (int ii = 1; ii < nx - 1; ii++)
      {
        for (int jj = 1; jj < ny - 1; jj++)
      {
        double rr = mesh[ii][jj] _MASS;
        double bx = mesh[ii][jj] _B_X;
        double by = mesh[ii][jj] _B_Y;
        double bz = mesh[ii][jj] _B_Z;

        if (isnan (fy[ii][jj] _ENER))
          {
            std::cout
          << __FUNCTION__ << " Flux (" << ii << "," << jj << ")";

            for (int hh = 0; hh < 8; ++hh)
          {
            std::cout << " fy " << fy[ii][jj].array[hh];
          }

            std::cout << std::endl;
            MPI_Abort (MPI_COMM_WORLD, 1222);
            exit (0);
          }

        // r and z are distances from (0,0) in global grid
        double myx = (double) myaddress[0] + ii;
        double myy = (double) myaddress[1] + jj;
        double r = myx * delta_x;
        double z = myy * delta_y;
        double GM = 1.0 / (eps * eps);
        //double etam=alpham * Va * H *exp(-2*zz*zz/(H*H));
        double etam = 0.1;
        double alpham = 0.1;
        // H is the thermal heightscale of the disk
        // cs is soundspeed at disk midplane
        // cs = sqrt (gamma *P/ rho);
        double cs = 0;
        cs = SoundSpeedMidplane[ii];
        // omegaK is keplerian velocity at disk midplane
        double omegaK = 0;
        omegaK = sqrt (GM / r * r * r);
        double H = eps * r;
        H = cs / omegaK;
        double Va =
          std::sqrt (bx * bx + by * by + bz * bz) / std::sqrt (rr);
        etam = alpham * H * Va * exp (-2 * z * z / (H * H));
        double chim = 1.0;
        double etamdash = chim * etam;
        double etaR = etamdash;
        double etaZ = etamdash;
        double etaPhi = etam;

        // Resistive Terms (See Toth 2000)
        // add a -B X (eta J)
        fx[ii][jj] _ENER +=
          (by * etaPhi * xjz[ii][jj] - bz * etaZ * xjy[ii][jj]);
        fy[ii][jj] _ENER +=
          (-bx * etaPhi * yjz[ii][jj] + bz * etaR * yjx[ii][jj]);
        if (isnan (fy[ii][jj] _ENER))
          {
            std::cout
          << __FUNCTION__
          << " Flux ("
          << ii << "," << jj
          << ")"
          << " bx " << bx
          << " by " << by
          << " bz " << bz
          << " rr " << rr
          << " etaPhi " << etaPhi
          << " yjz " << yjz[ii][jj]
          << " etaR " << etaR << " yjx " << yjx[ii][jj] << std::endl;
            MPI_Abort (MPI_COMM_WORLD, 1222);
            exit (0);
          }



      }
      }
#endif

    faceEx = 0.0;
    faceEy = 0.0;
    faceEx = 9.0E99;
    faceEy = 9.0e99;

    for (int ii = 2; ii < nx - 1; ii++) {
        for (int jj = 2; jj < ny; jj++) {
            faceEx[ii][jj] = fx[ii][jj]_B_Y;
        }
    }

    for (int ii = 2; ii < nx; ii++) {
        for (int jj = 2; jj < ny - 1; jj++) {
            faceEy[ii][jj] = fy[ii][jj]_B_X;
        }
    }

    return 0;
}
