/* $Id: parsecond.cpp,v 1.6 2006-11-16 13:48:07 gmurphy Exp $  */

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
parsecond(TNT::Array2D<unk> mesh,
          TNT::Array2D<flux> fx,
          TNT::Array2D<flux> fy,
          TNT::Array2D<flux> SecondOrdCorr,
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
          double delta_x, double delta_y, int *myaddress, MPI_Comm Cart_comm) {

    double pmin = 1e-5;
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



    // Sweep for fluxes
    for (int ii = 2; ii < nx - 1; ii++) {
        for (int jj = 2; jj < ny - 1; jj++) {

            double rr = mesh[ii][jj]_MASS;
            if (rr < 0) {
                std::cout
                        << __FUNCTION__ << ": "
                        << "Negative density at"
                        << " (" << ii << " , " << jj << " )" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1234);
                exit(0);
            }
            double px = mesh[ii][jj]_MOMX;
            double py = mesh[ii][jj]_MOMY;
            double pz = mesh[ii][jj]_MOMZ;
            double et = mesh[ii][jj]_ENER;
            double bx = mesh[ii][jj]_B_X;
            double by = mesh[ii][jj]_B_Y;
            double bz = mesh[ii][jj]_B_Z;

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


            if (isnan(p)) {
                std::cout << "Negative pressure " << std::endl;
                std::cout
                        << __FUNCTION__ << ": "
                        << " ("
                        << ii << "," << jj
                        << ") bx " << bx << " by " << by << " bz " << bz << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 44444);
                exit(0);
            }

            // Work out current at each cell face and add to energy flux
            // xflux current
            double idx = 1.0 / (delta_x);
            double idy = 1.0 / (delta_y);


            xjx[ii][jj] =
                    0.5 * idy * ((mesh[ii][jj + 1]_B_Z + mesh[ii - 1][jj + 1]_B_Z) -
                                 (mesh[ii][jj - 1]_B_Z + mesh[ii - 1][jj - 1]_B_Z));

            if (isnan(xjx[ii][jj])) {
                std::cout
                        << "xjx"
                        << " " << delta_y
                        << " " << idy
                        << " " << delta_x
                        << " " << idx
                        << " " << mesh[ii][jj + 1]_B_Z
                        << " " << mesh[ii - 1][jj + 1]_B_Z << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 666666);
                exit(0);
            }

            xjy[ii][jj] =
                    (-idx) * ((mesh[ii][jj]_B_Z - mesh[ii - 1][jj]_B_Z));
            xjz[ii][jj] =
                    (idx) * ((mesh[ii][jj]_B_Y - mesh[ii - 1][jj]_B_Y)) -
                    0.5 * idy * ((mesh[ii][jj + 1]_B_X + mesh[ii - 1][jj + 1]_B_X) -
                                 (mesh[ii][jj - 1]_B_X + mesh[ii - 1][jj - 1]_B_X));

            // yflux current
            yjx[ii][jj] = idy * (mesh[ii][jj]_B_Z - mesh[ii][jj - 1]_B_Z);
            yjy[ii][jj] =
                    0.5 * (-idx) *
                    ((mesh[ii + 1][jj]_B_Z + mesh[ii + 1][jj - 1]_B_Z) -
                     (mesh[ii - 1][jj]_B_Z + mesh[ii - 1][jj - 1]_B_Z));
            yjz[ii][jj] =
                    0.5 * (idx) *
                    ((mesh[ii + 1][jj]_B_Y + mesh[ii + 1][jj - 1]_B_Y) -
                     (mesh[ii - 1][jj]_B_Y + mesh[ii - 1][jj - 1]_B_Y)) -
                    idy * (mesh[ii][jj]_B_X - mesh[ii][jj - 1]_B_X);


            double lstate[8];
            double rstate[8];
            double iflux[8];
            double Res_state[8];
            int timestep = 22;
            double unused = 0;


            int rc = 0;


            int loc[2];
            double rz[2];
            loc[0] = (double) myaddress[0] + ii;
            loc[1] = (double) myaddress[1] + jj;
            double llprim[8];
            double leftprim[8];
            double rightprim[8];
            double rrprim[8];
            double llstate[8];
            double rrstate[8];


            rc = ctop(mesh[ii - 2][jj].array, llprim, rz);
            llprim[4] = std::max(llprim[4], pmin);
            errcheck(rc, myaddress, ii - 2, jj);

            rc = ctop(mesh[ii - 1][jj].array, leftprim, rz);
            leftprim[4] = std::max(leftprim[4], pmin);
            errcheck(rc, myaddress, ii - 1, jj);

            rc = ctop(mesh[ii][jj].array, rightprim, rz);
            rightprim[4] = std::max(rightprim[4], pmin);
            errcheck(rc, myaddress, ii, jj);

            rc = ctop(mesh[ii + 1][jj].array, rrprim, rz);
            rrprim[4] = std::max(rrprim[4], pmin);
            errcheck(rc, myaddress, ii + 1, jj);


            // Limiter
            //


            double myxi = myaddress[0] + (ii - NGC);
            double rgi =
                    (myxi * myxi + myxi + 1.0 / 3.0) * delta_x / (myxi + 0.5);
            double myxp = myaddress[0] + (ii - NGC) + 1;
            double rgii =
                    (myxp * myxp + myxp + 1.0 / 3.0) * delta_x / (myxp + 0.5);


            for (int hh = 0; hh < 8; hh++) {
                double left = llprim[hh];
                double mid = leftprim[hh];
                double right = rightprim[hh];
                double slope1 = LIMITER((mid - left), (right - mid));
                //
                left = leftprim[hh];
                mid = rightprim[hh];
                right = rrprim[hh];
                double slope2 = LIMITER((mid - left), (right - mid));
                //
#ifdef CYLINDRICAL
                //leftprim[hh] = leftprim[hh] + ((myxi + 1) * delta_x - rgi) * slope1;
                //rightprim[hh] = rightprim[hh] - (rgii - (myxi + 1) * delta_x) * slope2;
                leftprim[hh] = leftprim[hh] + 0.5 * delta_x * slope1;
                rightprim[hh] = rightprim[hh] - 0.5 * delta_x * slope2;
                SecondOrdCorr[ii][jj].array[hh] = slope1;
#else
                leftprim[hh] = leftprim[hh] + 0.5 * delta_x * slope1;
                rightprim[hh] = rightprim[hh] - 0.5 * delta_x * slope2;
#endif
            }


/*
 			hh=4;
	    {
	      double left = llprim[hh];
	      double mid = leftprim[hh];
	      double right = rightprim[hh];
	      double slope1 = minmod ((mid - left), (right - mid));
	      //
	      left = leftprim[hh];
	      mid = rightprim[hh];
	      right = rrprim[hh];
	      double slope2 = minmod ((mid - left), (right - mid));
	      //
#ifdef CYLINDRICAL
	      leftprim[hh] =
		leftprim[hh] + ((myxi + 1) * delta_x - rgi) * slope1;
	      rightprim[hh] =
		rightprim[hh] - (rgii - (myxi + 1) * delta_x) * slope2;
	      SecondOrdCorr[ii][jj].array[hh] = slope1;
#else
	      leftprim[hh] = leftprim[hh] + 0.5 * delta_x * slope1;
	      rightprim[hh] = rightprim[hh] - 0.5 * delta_x * slope2;
#endif
	    }
		 */

            leftprim[4] = std::max(leftprim[4], pmin);
            rightprim[4] = std::max(rightprim[4], pmin);


#ifdef HLLD
#ifdef STAGGER_MESH
            leftprim[5] = faceBx[ii][jj];
            rightprim[5] = faceBx[ii][jj];
            /*
            if ( fabs(faceBx[ii][jj]) > 1e10)
            {
                std::cout
                << " ( " << ii << " , " << jj << " ) "
                << " faceBx " << faceBx[ii][jj]
                << std::endl;
                }
                */
#endif
            rc = hlld(leftprim, rightprim, fx[ii][jj].array, Res_state, timestep, &unused, 1, loc);
#ifdef LAXFRIEDRICHS
            rc = lf (leftprim, rightprim, fx[ii][jj].array, Res_state, timestep, &unused, 1, loc);
#endif
            errcheck(rc, myaddress, ii, jj);
#endif

#ifdef ROE
            rc = ptoc (leftprim, lstate);
            rc = ptoc (rightprim, rstate);
            rc =
              riemann (lstate, rstate,
                   fx[ii][jj].array, Res_state, timestep, &unused, 1);
#endif

            llstate[0] = mesh[ii][jj - 2].array[0];
            llstate[1] = mesh[ii][jj - 2].array[2];
            llstate[2] = mesh[ii][jj - 2].array[3];
            llstate[3] = mesh[ii][jj - 2].array[1];
            llstate[4] = mesh[ii][jj - 2].array[4];
            llstate[5] = mesh[ii][jj - 2].array[6];
            llstate[6] = mesh[ii][jj - 2].array[7];
            llstate[7] = mesh[ii][jj - 2].array[5];

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


            rrstate[0] = mesh[ii][jj + 1].array[0];
            rrstate[1] = mesh[ii][jj + 1].array[2];
            rrstate[2] = mesh[ii][jj + 1].array[3];
            rrstate[3] = mesh[ii][jj + 1].array[1];
            rrstate[4] = mesh[ii][jj + 1].array[4];
            rrstate[5] = mesh[ii][jj + 1].array[6];
            rrstate[6] = mesh[ii][jj + 1].array[7];
            rrstate[7] = mesh[ii][jj + 1].array[5];

            rc = ctop(llstate, llprim, rz);
            llprim[4] = std::max(llprim[4], pmin);
            errcheck(rc, myaddress, ii, jj - 1);
            rc = ctop(lstate, leftprim, rz);
            leftprim[4] = std::max(leftprim[4], pmin);
            errcheck(rc, myaddress, ii, jj - 1);
            rc = ctop(rstate, rightprim, rz);
            rightprim[4] = std::max(rightprim[4], pmin);
            errcheck(rc, myaddress, ii, jj);
            rc = ctop(rrstate, rrprim, rz);
            rrprim[4] = std::max(rrprim[4], pmin);
            errcheck(rc, myaddress, ii, jj);


            for (int hh = 0; hh < 8; hh++) {
                double left = llprim[hh];
                double mid = leftprim[hh];
                double right = rightprim[hh];
                double slope1 = LIMITER((mid - left), (right - mid));
                //
                left = leftprim[hh];
                mid = rightprim[hh];
                right = rrprim[hh];
                double slope2 = LIMITER((mid - left), (right - mid));
                //
                leftprim[hh] = leftprim[hh] + 0.5 * delta_y * slope1;
                rightprim[hh] = rightprim[hh] - 0.5 * delta_y * slope2;
            }

            leftprim[4] = std::max(leftprim[4], pmin);
            rightprim[4] = std::max(rightprim[4], pmin);
#ifdef HLLD
#ifdef STAGGER_MESH
            leftprim[5] = faceBy[ii][jj];
            rightprim[5] = faceBy[ii][jj];
            /*
            if ( fabs(faceBy[ii][jj]) > 1e10)
            {
                std::cout
                << " ( " << ii << " , " << jj << " ) "
                << " " << faceBy[ii][jj]
                << std::endl;
                }
                */

#endif
            rc = hlld(leftprim, rightprim, iflux, Res_state, timestep, &unused, 2, loc);
#ifdef LAXFRIEDRICHS
            rc = lf (leftprim, rightprim, iflux, Res_state, timestep, &unused, 1, loc);
#endif

#endif

#ifdef ROE
            rc = ptoc (leftprim, lstate);
            rc = ptoc (rightprim, rstate);
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

            if (isnan(fy[ii][jj]_ENER)) {
                std::cout
                        << __FUNCTION__ << " Flux (" << ii << "," << jj << ")";

                for (int hh = 0; hh < 8; ++hh) {
                    std::cout << " fy " << fy[ii][jj].array[hh];
                }

                std::cout << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1222);
                exit(0);
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
            omegaK = sqrt(GM / r * r * r);
            double H = eps * r;
            H = cs / omegaK;
            double Va =
                    std::sqrt(bx * bx + by * by + bz * bz) / std::sqrt(rr);
            etam = alpham * H * Va * exp(-2 * z * z / (H * H));
            double chim = 1.0;
            double etamdash = chim * etam;
            double etaR = etamdash;
            double etaZ = etamdash;
            double etaPhi = etam;

            // Resistive Terms
            // Toth 2000 J. Comp. Phys.
            // add a -B X (eta J)
#ifdef RESISTIVE
            fx[ii][jj] _ENER +=
              (by * etaPhi * xjz[ii][jj] - bz * etaZ * xjy[ii][jj]);
            fy[ii][jj] _ENER +=
              (-bx * etaPhi * yjz[ii][jj] + bz * etaR * yjx[ii][jj]);
#endif
            if (isnan(fy[ii][jj]_ENER)) {
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
                MPI_Abort(MPI_COMM_WORLD, 1222);
                exit(0);
            }


        }
    }

    faceEx = 9e99;
    faceEy = 9e99;
    faceEx = 0.0;
    faceEy = 0.0;

    for (int ii = 1; ii < nx - 1; ii++) {
        for (int jj = 1; jj < ny; jj++) {
            faceEx[ii][jj] = fx[ii][jj]_B_Y;
        }
    }

    for (int ii = 1; ii < nx; ii++) {
        for (int jj = 1; jj < ny - 1; jj++) {
            faceEy[ii][jj] = fy[ii][jj]_B_X;
        }
    }

}

