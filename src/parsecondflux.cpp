/* $Id: parsecondflux.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */

#include "mpi.h"
#include "tnt.h"
#include "out.h"
#include "PhysConsts.h"
#include "parflux.h"
int
parflux (TNT::Array2D < unk > mesh,
	 TNT::Array2D < flux > fx,
	 TNT::Array2D < flux > fy,
	 TNT::Array2D < double >xjx,
	 TNT::Array2D < double >xjy,
	 TNT::Array2D < double >xjz,
	 TNT::Array2D < double >yjx,
	 TNT::Array2D < double >yjy,
	 TNT::Array2D < double >yjz,
	 double *delta_x, double *delta_y, int *myaddress, MPI_Comm Cart_comm)
{

  double eps = 0.1;
  double gammam1 = PhysConsts::gamma - 1.;
  double gammam1i = 1. / gammam1;

  int nx = mesh.dim1 ();
  int ny = mesh.dim2 ();
  // Current arrays

  int myid = 99999;
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  int ndims = 2;
  int coords[2];


  int remain_dims[] = { 0, 1 };
  MPI_Comm Sub_comm;
  MPI_Cart_sub (Cart_comm, remain_dims, &Sub_comm);
  int lowest = 99999;
  coords[1] = 1;
  MPI_Cart_coords (Cart_comm, myid, ndims, coords);
  MPI_Cart_rank (Cart_comm, coords, &lowest);
  // Check my position in the Cartesian grid
  // Then bcast my id to the Sub_comm
  // Then all print out the id
  // Grab the coords of the processor
  double tmpbuf;
  tmpbuf = myid;
  TNT::Array1D < double >pressure_buf (nx);
  MPI_Bcast (&tmpbuf, 1, MPI_DOUBLE, 0, Sub_comm);

#ifdef DEBUG_SUB_COMM
  std::cout
    << "myid = " << myid
    << "coords = (" << coords[0] << "," << coords[1] << ")"
    << " lowest = " << tmpbuf << std::endl;
#endif



  // Sweep for fluxes 
  for (int ii = 1; ii <= nx - 2; ii++)
    {
      for (int jj = 1; jj <= ny - 2; jj++)
	{

	  double rr = mesh[ii][jj] _MASS;
	  if (rr < 0)
	  {
		  std::cout 
			  << __FUNCTION__ << ": "
			  << "Negative density at"
			  << " (" << ii << " , " << jj <<" )"
			  <<std::endl;
		  MPI_Abort(MPI_COMM_WORLD, 1234);
		  exit(0);
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
	  double ke = 0.5 * rr * (vx * vx + vy * vy);
	  double b2 = 0.5 * (bx * bx + by * by + bz * bz);
	  double vdotb = (vx * bx + vy * by + vz * bz);
	  double p = et - ke - b2;
	  double ptot = et - ke;
	  p = p * gammam1;


	  if (isnan (p))
	    {
	      std::cout
		<< __FUNCTION__
		<< " ("
		<< ii << "," << jj
		<< ") bx " << bx << " by " << by << " bz " << bz << std::endl;
			MPI_Abort (MPI_COMM_WORLD, 44444);
			exit(0);
	    }

	  // Work out current at each cell face and add to energy flux
	  // xflux current
	  double idx = 1.0 / (*delta_x);
	  double idy = 1.0 / (*delta_y);


	  xjx[ii][jj] =
	    0.5 * idy * ((mesh[ii][jj + 1] _B_Z + mesh[ii - 1][jj + 1] _B_Z) -
			 (mesh[ii][jj - 1] _B_Z + mesh[ii - 1][jj - 1] _B_Z));

	  if (isnan (xjx[ii][jj]))
	    {
	      std::cout
		<< "xjx"
		<< " " << *delta_y
		<< " " << idy
		<< " " << *delta_x
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



	  fx[ii][jj] _MASS = rr * vx;
	  fx[ii][jj] _MOMX = rr * vx * vx + ptot - bx * bx;
	  fx[ii][jj] _MOMY = rr * vx * vy - bx * by;
	  fx[ii][jj] _MOMZ = rr * vx * vz - bx * bz;
	  fx[ii][jj] _ENER = (et + ptot) * vx - bx * vdotb;
	  fx[ii][jj] _B_X = 0;
	  fx[ii][jj] _B_Y = by * bx - bx * vy;
	  fx[ii][jj] _B_Z = bz * vx - bx * vz;

	  fy[ii][jj] _MASS = rr * vy;
	  fy[ii][jj] _MOMX = rr * vx * vy - by * bx;
	  fy[ii][jj] _MOMY = rr * vy * vy + ptot - by * by;
	  fy[ii][jj] _MOMZ = rr * vz * vy - by * bz;
	  fy[ii][jj] _ENER = (et + ptot) * vy - by * vdotb;
	  fy[ii][jj] _B_X = bx * vy - by * vx;
	  fy[ii][jj] _B_Y = 0;
	  fy[ii][jj] _B_Z = bz * vy - by * vz;


	  if (isnan (fy[ii][jj] _ENER))
	    {
	      std::cout
		<< __FUNCTION__
		<< " Flux ("
		<< ii << "," << jj
		<< ") et " << et
		<< " ptot " << ptot
		<< " p " << p
		<< " vy " << vy
		<< " by " << by
		<< " vdotb " << vdotb << " rr " << rr << std::endl;
	      MPI_Abort (MPI_COMM_WORLD, 1222);
	      exit (0);
	    }

	  // Gravity Term
	  // r and z are distances from (0,0) in global grid
	  double myx = (double) myaddress[0] + ii;
	  double myy = (double) myaddress[1] + jj;
	  double r = myx * *delta_x;
	  double z = myy * *delta_y;
#define GRAVITY
#ifdef GRAVITY
	  double GM = 1.0 / (eps * eps);
	  double gravity =
	    (-GM) / std::sqrt (std::pow (r, 2.0) + std::pow (z, 2.0));
	  double grav_en = mesh[ii][jj] _MOMX * gravity;
	  //      std::cout << " r " << r << " z " << z << std::endl;
	  gravity = mesh[ii][jj] _MASS *gravity;
	  fx[ii][jj] _MOMX += gravity;
	  fx[ii][jj] _ENER += grav_en;

	  grav_en = mesh[ii][jj] _MOMY *gravity;
	  gravity = mesh[ii][jj] _MASS *gravity;
	  fy[ii][jj] _MOMY += gravity;
	  fy[ii][jj] _ENER += grav_en;
#else
	  std::cout << "Gravity disabled" << std::endl;

#endif
	  //double etam=alpham * Va * H *exp(-2*zz*zz/(H*H));
	  double etam = 0.1;
	  double alpham = 0.1;
	  // H is the thermal heightscale of the disk
	  // cs is soundspeed at disk midplane
	  // cs = sqrt (gamma *P/ rho);
	  // double cs=0;
	  // omegaK is keplerian velocity at disk midplane
	  double omegaK = 0;
	  omegaK = sqrt (GM / r * r * r);
	  double H = eps * r;
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

}
