/* $Id: parboundary.cpp,v 1.5 2006-11-16 13:48:07 gmurphy Exp $  */

#include "problem.h"
#include "physics.h"
#include "parboundary.h"
int
parboundary (
TNT::Array2D < unk > mesh,
TNT::Array2D < unk > recvmesh,
	     TNT::Array2D < double >faceBx,
	     TNT::Array2D < double >faceBy,
	     int myNorth, int mySouth, int myEast, int myWest, int myid,
	     MPI_Comm Cart_comm)
{


  MPI_Status status;

  int nx = mesh.dim1 ();
  int ny = mesh.dim2 ();
  int fnx = faceBx.dim1 ();
  int fny = faceBy.dim2 ();
  int ne = NE;

  unk test;


  double gammam1 = PhysConsts::gamma - 1;
  double gammam1i = 1.0 / gammam1;


  int myaddress[] = { 0, 0 };

  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  int ndims = 2;
  int coords[2];
  MPI_Cartdim_get (Cart_comm, &ndims);
  MPI_Cart_coords (Cart_comm, myid, ndims, coords);
  myaddress[0] = (nx - 2 * NGC) * coords[0];
  myaddress[1] = (ny - 2 * NGC) * coords[1];


  MPI_Datatype MPI_unk_type;
  MPI_Datatype MpiColumnType;
  MPI_Datatype MpiRowType;
  MPI_Datatype MpiFaceColumnType;
  MPI_Datatype MpiFaceRowType;
  MPI_Datatype MpiyFaceColumnType;
  MPI_Datatype MpiyFaceRowType;



  MPI_Aint disps[3] = { 0, (unsigned long int) test.array - (unsigned long int) &test.temperature, sizeof (test) };	/* guessing... */
  int blks[3] = { 1, ne, 1 };
  MPI_Datatype types[3] = { MPI_INT, MPI_DOUBLE, MPI_UB };


  MPI_Type_struct (3, blks, disps, types, &MPI_unk_type);
  MPI_Type_commit (&MPI_unk_type);


  MPI_Type_vector (nx,		/* # column elements */
		   1,		/* 1 column only */
		   ny,		/* skip ny elements */
		   MPI_unk_type,	/* elements are double */
		   &MpiColumnType);	/* MPI derived datatype */

  MPI_Type_vector (1,		/* # row elements */
		   ny,		/* 1 row only */
		   0,		/* skip ny elements */
		   //    myrow,   /* elements are double */
		   MPI_unk_type,	/* elements are double */
		   &MpiRowType);	/* MPI derived datatype */


  MPI_Type_vector (nx + 1,	/* # column elements */
		   1,		/* 1 column only */
		   ny + 1,	/* skip ny elements */
		   MPI_DOUBLE,	/* elements are double */
		   &MpiFaceColumnType);	/* MPI derived datatype */

  MPI_Type_vector (1,		/* # row elements */
		   ny + 1,	/* 1 row only */
		   0,		/* skip ny elements */
		   //    myrow,   /* elements are double */
		   MPI_DOUBLE,	/* elements are double */
		   &MpiFaceRowType);	/* MPI derived datatype */



  MPI_Type_vector (nx + 1,	/* # column elements */
		   1,		/* 1 column only */
		   ny + 1,	/* skip ny elements */
		   MPI_DOUBLE,	/* elements are double */
		   &MpiyFaceColumnType);	/* MPI derived datatype */

  MPI_Type_vector (1,		/* # row elements */
		   ny + 1,	/* 1 row only */
		   0,		/* skip ny elements */
		   //    myrow,   /* elements are double */
		   MPI_DOUBLE,	/* elements are double */
		   &MpiyFaceRowType);	/* MPI derived datatype */


  MPI_Type_commit (&MpiColumnType);
  MPI_Type_commit (&MpiRowType);
  MPI_Type_commit (&MpiFaceColumnType);
  MPI_Type_commit (&MpiFaceRowType);
  MPI_Type_commit (&MpiyFaceColumnType);
  MPI_Type_commit (&MpiyFaceRowType);



#define PARALLEL

#ifdef PARALLEL

  // Send Columns
  /* Sends to proc above */// Send him my second col from top
  MPI_Send (&mesh[0][ny - 4], 1, MpiColumnType, myNorth, myid,
	    MPI_COMM_WORLD);
  /* Sends to proc below */// Send him my second lowest col
  MPI_Send (&mesh[0][3], 1, MpiColumnType, mySouth, myid, MPI_COMM_WORLD);	/* Receives from proc below */
  // Receive into lowest col
  MPI_Recv (&mesh[0][0], 1, MpiColumnType, mySouth, mySouth, MPI_COMM_WORLD,
	    &status);
  /* Receives from proc above */// Receive into highest col
  MPI_Recv (&mesh[0][ny - 1], 1, MpiColumnType, myNorth, myNorth,
	    MPI_COMM_WORLD, &status);
  MPI_Barrier (MPI_COMM_WORLD);
  /* Sends to proc above */// Send him my second col from top
  MPI_Send (&mesh[0][ny - 3], 1, MpiColumnType, myNorth, myid,
	    MPI_COMM_WORLD);
  /* Sends to proc below */// Send him my second lowest col
  MPI_Send (&mesh[0][2], 1, MpiColumnType, mySouth, myid, MPI_COMM_WORLD);
  /* Receives from proc below */// Receive into lowest col
  MPI_Recv (&mesh[0][1], 1, MpiColumnType, mySouth, mySouth, MPI_COMM_WORLD,
	    &status);
  /* Receives from proc above */// Receive into highest col
  MPI_Recv (&mesh[0][ny - 2], 1, MpiColumnType, myNorth, myNorth,
	    MPI_COMM_WORLD, &status);
  MPI_Barrier (MPI_COMM_WORLD);



  // Send Rows
  /* Sends to proc W */
  MPI_Send (&mesh[3][0], 1, MpiRowType, myWest, myid, MPI_COMM_WORLD);
  /* Sends to proc E */
  MPI_Send (&mesh[nx - 4][0], 1, MpiRowType, myEast, myid, MPI_COMM_WORLD);
  /* Receives from proc E */
  MPI_Recv (&mesh[nx - 1][0], 1, MpiRowType, myEast, myEast, MPI_COMM_WORLD,
	    &status);
  /* Receives from proc W */
  MPI_Recv (&mesh[0][0], 1, MpiRowType, myWest, myWest, MPI_COMM_WORLD,
	    &status);
  MPI_Barrier (MPI_COMM_WORLD);
  // Send Rows
  /* Sends to proc W */
  MPI_Send (&mesh[2][0], 1, MpiRowType, myWest, myid, MPI_COMM_WORLD);
  /* Sends to proc E */
  MPI_Send (&mesh[nx - 3][0], 1, MpiRowType, myEast, myid, MPI_COMM_WORLD);
  /* Receives from proc E */
  MPI_Recv (&mesh[nx - 2][0], 1, MpiRowType, myEast, myEast, MPI_COMM_WORLD,
	    &status);
  /* Receives from proc W */
  MPI_Recv (&mesh[1][0], 1, MpiRowType, myWest, myWest, MPI_COMM_WORLD,
	    &status);
  MPI_Barrier (MPI_COMM_WORLD);

#else 


  for (int jj = 0; jj < ny; jj++)
  {
      mesh[0][jj]= mesh[nx-4][jj];
      mesh[1][jj]= mesh[nx-3][jj];
      mesh[nx-2][jj]= mesh[2][jj];
      mesh[nx-1][jj]= mesh[3][jj];
  }
    for (int ii = 0; ii < nx; ii++)
    {
      mesh[ii][0]= mesh[ii][ny-4];
      mesh[ii][1]= mesh[ii][ny-3];
      mesh[ii][ny-2]= mesh[ii][2];
      mesh[ii][ny-1]= mesh[ii][3];
       }

		
#endif


#ifdef EXCHANGE_FACE

  //========== FACE VARIABLES ===============

  // Send Columns
  // Sends to proc above // Send him my 4th col from top
  MPI_Send (&faceBx[0][fny - 4], 1, MpiFaceColumnType, myNorth, myid,
	    MPI_COMM_WORLD);
  // Sends to proc below // Send him my 4th lowest col
  MPI_Send (&faceBx[0][3], 1, MpiFaceColumnType, mySouth, myid,
	    MPI_COMM_WORLD);
  // Receives from proc below  // Receive into lowest col
  MPI_Recv (&faceBx[0][0], 1, MpiFaceColumnType, mySouth, mySouth,
	    MPI_COMM_WORLD, &status);
  // Receives from proc above // Receive into highest col
  MPI_Recv (&faceBx[0][fny - 1], 1, MpiFaceColumnType, myNorth, myNorth,
	    MPI_COMM_WORLD, &status);
  MPI_Barrier (MPI_COMM_WORLD);
  // Sends to proc above // Send him my 3rd highest col 
  MPI_Send (&faceBx[0][fny - 3], 1, MpiFaceColumnType, myNorth, myid,
	    MPI_COMM_WORLD);
  // Sends to proc below // Send him my 3rd lowest col
  MPI_Send (&faceBx[0][2], 1, MpiFaceColumnType, mySouth, myid,
	    MPI_COMM_WORLD);
  // Receives from proc below // Receive into 2nd lowest col
  MPI_Recv (&faceBx[0][1], 1, MpiFaceColumnType, mySouth, mySouth,
	    MPI_COMM_WORLD, &status);
  // Receives from proc above // Receive into 2nd highest col
  MPI_Recv (&faceBx[0][fny - 2], 1, MpiFaceColumnType, myNorth, myNorth,
	    MPI_COMM_WORLD, &status);
  MPI_Barrier (MPI_COMM_WORLD);

  // Send Rows
  // Sends to proc W 
  MPI_Send (&faceBx[3][0], 1, MpiFaceRowType, myWest, myid, MPI_COMM_WORLD);
  // Sends to proc E 
  MPI_Send (&faceBx[fnx - 4][0], 1, MpiFaceRowType, myEast, myid,
	    MPI_COMM_WORLD);
  // Receives from proc E 
  MPI_Recv (&faceBx[fnx - 1][0], 1, MpiFaceRowType, myEast, myEast,
	    MPI_COMM_WORLD, &status);
  // Receives from proc W 
  MPI_Recv (&faceBx[0][0], 1, MpiFaceRowType, myWest, myWest, MPI_COMM_WORLD,
	    &status);
  MPI_Barrier (MPI_COMM_WORLD);
  // Send FaceRows
  // Sends to proc W 
  MPI_Send (&faceBx[2][0], 1, MpiFaceRowType, myWest, myid, MPI_COMM_WORLD);
  // Sends to proc E 
  MPI_Send (&faceBx[fnx - 3][0], 1, MpiFaceRowType, myEast, myid,
	    MPI_COMM_WORLD);
  // Receives from proc E 
  MPI_Recv (&faceBx[fnx - 2][0], 1, MpiFaceRowType, myEast, myEast,
	    MPI_COMM_WORLD, &status);
  // Receives from proc W 
  MPI_Recv (&faceBx[1][0], 1, MpiFaceRowType, myWest, myWest, MPI_COMM_WORLD,
	    &status);
  MPI_Barrier (MPI_COMM_WORLD);


  // Send Columns
  // Sends to proc above // Send him my second col from top
  MPI_Send (&faceBy[0][fny - 4], 1, MpiFaceColumnType, myNorth, myid,
	    MPI_COMM_WORLD);
  // Sends to proc below // Send him my second lowest col
  MPI_Send (&faceBy[0][3], 1, MpiFaceColumnType, mySouth, myid, MPI_COMM_WORLD);	// Receives from proc below 
  // Receive into lowest col
  MPI_Recv (&faceBy[0][0], 1, MpiFaceColumnType, mySouth, mySouth,
	    MPI_COMM_WORLD, &status);
  // Receives from proc above // Receive into highest col
  MPI_Recv (&faceBy[0][fny - 1], 1, MpiFaceColumnType, myNorth, myNorth,
	    MPI_COMM_WORLD, &status);
  MPI_Barrier (MPI_COMM_WORLD);
  // Sends to proc above // Send him my second col from top
  MPI_Send (&faceBy[0][fny - 3], 1, MpiFaceColumnType, myNorth, myid,
	    MPI_COMM_WORLD);
  // Sends to proc below // Send him my second lowest col
  MPI_Send (&faceBy[0][2], 1, MpiFaceColumnType, mySouth, myid,
	    MPI_COMM_WORLD);
  // Receives from proc below // Receive into lowest col
  MPI_Recv (&faceBy[0][1], 1, MpiFaceColumnType, mySouth, mySouth,
	    MPI_COMM_WORLD, &status);
  // Receives from proc above // Receive into highest col
  MPI_Recv (&faceBy[0][fny - 2], 1, MpiFaceColumnType, myNorth, myNorth,
	    MPI_COMM_WORLD, &status);
  MPI_Barrier (MPI_COMM_WORLD);

  // Send Rows
  // Sends to proc W 
  MPI_Send (&faceBy[3][0], 1, MpiFaceRowType, myWest, myid, MPI_COMM_WORLD);
  // Sends to proc E 
  MPI_Send (&faceBy[fnx - 4][0], 1, MpiFaceRowType, myEast, myid,
	    MPI_COMM_WORLD);
  // Receives from proc E 
  MPI_Recv (&faceBy[fnx - 1][0], 1, MpiFaceRowType, myEast, myEast,
	    MPI_COMM_WORLD, &status);
  // Receives from proc W 
  MPI_Recv (&faceBy[0][0], 1, MpiFaceRowType, myWest, myWest, MPI_COMM_WORLD,
	    &status);
  MPI_Barrier (MPI_COMM_WORLD);
  // Send FaceRows
  // Sends to proc W 
  MPI_Send (&faceBy[2][0], 1, MpiFaceRowType, myWest, myid, MPI_COMM_WORLD);
  // Sends to proc E 
  MPI_Send (&faceBy[fnx - 3][0], 1, MpiFaceRowType, myEast, myid,
	    MPI_COMM_WORLD);
  // Receives from proc E 
  MPI_Recv (&faceBy[fnx - 2][0], 1, MpiFaceRowType, myEast, myEast,
	    MPI_COMM_WORLD, &status);
  // Receives from proc W 
  MPI_Recv (&faceBy[1][0], 1, MpiFaceRowType, myWest, myWest, MPI_COMM_WORLD,
	    &status);
  MPI_Barrier (MPI_COMM_WORLD);
#endif

#ifdef MAES
  if (myNorth == MPI_PROC_NULL)
    {
      // Upper y boundary is freeflow
      for (int i = 0; i < nx; ++i)
	{
	  recvmesh[i][ny - 2] = mesh[i][ny - 3];
	  recvmesh[i][ny - 1] = mesh[i][ny - 2];
	}
    }



		

  if (mySouth == MPI_PROC_NULL)
    {
      //  Lower yboundary is reflecting
      for (int i = 0; i < nx; ++i)
	{
	  recvmesh[i][1] = mesh[i][2];
	  recvmesh[i][0] = mesh[i][3];
	  recvmesh[i][0] _MOMY = (-mesh[i][0] _MOMY);
	  recvmesh[i][1] _MOMY = (-mesh[i][1] _MOMY);
	}
    }

  if (myEast == MPI_PROC_NULL)
    {
      // Upper xboundary is freeflow
      for (int i = 0; i < nx; ++i)
	{
	  recvmesh[nx - 2][i] = mesh[nx - 3][i];
	  recvmesh[nx - 1][i] = mesh[nx - 2][i];
	}
      // Except below the disk height where material is injected
      if (mySouth == MPI_PROC_NULL)
	{
	  for (int i = 0; i < 40; ++i)
	    {
	      double r = 40.0;
	      double r2 = r * r;
	      double z = i;
	      double r0 = 4.0;
	      double r02 = r0 * r0;
	      double eps = 0.1;
	      double z2 = z * z;
	      double h = eps * r;
			if (z < h)
			{
	      double h2 = h * h;
	      double rho =
		std::max (1e-6,
			  (pow (4, 1.5) / pow ((4 * 4 + 40 * 40), 0.75))) *
		std::pow (std::
			  max (1e-6, (1 - 0.5 * (gammam1) * i / 0.1 * 40)),
			  gammam1i);


	      recvmesh[nx - 2][i] _MASS = rho;

	      // Rotation Profile
	      double vtheta =
		(1 -
		 eps * eps) * sqrt (r0) * exp (-2 * z2 / h2) / (eps *
								pow ((r02 +
								      r2),
								     0.25));
	      recvmesh[nx - 2][i] _MOMZ = mesh[nx - 2][i] _MASS *vtheta;



	      // ms is a parameter smaller than unity 
	      // which then ensures an initial subsonic poloidal 
	      // inflow
	      double ms = 0.3;


	      // Poloidal Velocity
	      double vr =
		(-ms) * sqrt (r0) * exp (-2 * z2 / h2) /
		(pow ((r02 + r2), 0.25));
	      recvmesh[nx - 2][i] _MOMX = mesh[nx - 2][i] _MASS *vr;


	      double vz = vr * z / r;
	      recvmesh[nx - 2][i] _MOMY = mesh[nx - 2][i] _MOMX *z / r;


	      // plasma beta parameter measuring the ratio of the 
	      // thermal pressure to the magnetic pressure at Z = 0
	      double beta = 1.0;

	      // Magnetic Field
	      double bz =
		pow (r0, 2.5) / (pow ((r02 + r2), 1.25) * sqrt (beta));
	      bz = 0.5 * (faceBy[nx - 2][i] + faceBy[nx - 2][i + 1]);
	      bz=mesh[nx - 3][i] _B_Y ;
	      double br = 0;
	      double btheta = 0;


	      // Pressure = rho ^ (5/3)
	      double pressure = pow (rho, PhysConsts::gamma);



	      recvmesh[nx - 2][i] _B_X = 0.;
	      recvmesh[nx - 2][i] _B_Y = bz;
	      recvmesh[nx - 2][i] _B_Z = 0.;
	      double vv2 = vz * vz + vr * vr + vtheta * vtheta;
	      double b2 = bz * bz + br * br + btheta * btheta;
	      recvmesh[nx - 2][i] _ENER =
		(0.5 * rho * vv2) + (0.5 * b2) + pressure * gammam1i;

	      recvmesh[nx - 1][i] = mesh[nx - 2][i];
			}

	    }
	}


    }

  if (myWest == MPI_PROC_NULL)
    {
      // Lower xboundary is reflecting axisymmetric
      for (int i = 0; i < nx - 0; ++i)
	{
	  recvmesh[1][i] = mesh[2][i];
	  recvmesh[0][i] = mesh[3][i];
//#define DEBUG_SYM
#ifndef DEBUG_SYM
	  recvmesh[0][i] _MOMX = (-mesh[0][i] _MOMX);
	  recvmesh[1][i] _MOMX = (-mesh[1][i] _MOMX);
	  recvmesh[0][i] _MOMZ = (-mesh[0][i] _MOMZ);
	  recvmesh[1][i] _MOMZ = (-mesh[1][i] _MOMZ);
	  recvmesh[0][i] _B_X = (-mesh[0][i] _B_X);
	  recvmesh[1][i] _B_X = (-mesh[1][i] _B_X);
	  recvmesh[0][i] _B_Z = (-mesh[0][i] _B_Z);
	  recvmesh[1][i] _B_Z = (-mesh[1][i] _B_Z);
#else
	  //std::cout << "Debugging symmetric boundary" << std::endl;
#endif
	}
    }


  // Internal Boundary
  // No material from internal boundary may enter the grid
  if (mySouth == MPI_PROC_NULL)
    {
      if (myWest == MPI_PROC_NULL)
	{
		int diskheight=8;
		int boxwidth=8;
	  for (int j = 0; j < diskheight; ++j)
	    {
	      for (int i = 0; i < boxwidth; ++i)
		{
		  recvmesh[i][diskheight-1] = mesh[i][diskheight];
		  recvmesh[i][diskheight-2] = mesh[i][diskheight];
		  recvmesh[i][j] = mesh[i][diskheight];
		  //recvmesh[i][j] _MOMX = std::min (0.0, mesh[i][j] _MOMX);
		  //recvmesh[i][j] _MOMY = std::min (0.0, mesh[i][j] _MOMY);
		}
	    }
	  for (int j = 0; j < diskheight; ++j)
	    {
		  	recvmesh[boxwidth-1][j] = mesh[boxwidth][j];
		  	recvmesh[boxwidth-2][j] = mesh[boxwidth][j];
		  //recvmesh[boxwidth-1][j] _MOMX = std::min (0.0, mesh[boxwidth-1][j] _MOMX);
		  //recvmesh[boxwidth-2][j] _MOMY = std::min (0.0, mesh[boxwidth-2][j] _MOMY);
		 }
	}
    }
#endif
#ifdef JET


  if (mySouth == MPI_PROC_NULL && myWest == MPI_PROC_NULL)
    {
      for (int i = 2; i < nx; ++i)
	{
	  recvmesh[i][1] = mesh[i][2];
	  recvmesh[i][0] = mesh[i][3];
	  mesh[i][0] _MOMY = (-mesh[i][0] _MOMY);
	  mesh[i][1] _MOMY = (-mesh[i][1] _MOMY);

	  double rho = 0.25;
	  double vx = 0.0;
	  double vy = 0.0;


	  double myx = (i - NGC);

	      double vel = 10.0;
//#define CHECK
#ifndef CHECK
	  if (myx < 40)
	    {
	      double vel = 10.0;
	      vy = vel;
	      if (myx > 30)
		{
		  vy = (40 - 30) * (myx - 40) / (0 - vel);
		}
	      // vy = cos (2*PhysConsts::pi *myx /40.0)
	    }
#else
	      vy = vel;
#endif



	  double vz = 0.0;
	  double bx = mesh[i][2] _B_X;
	  double by = mesh[i][2] _B_Y;
	  double bz = mesh[i][2] _B_Z;
	  double pressure = 0.0004;

	  const double gammam1i = 1 / (PhysConsts::gamma - 1);

	  double bsqr = 0.5 * (bx * bx + by * by + bz * bz);
	  double ke = 0.5 * rho * (vx * vx + vy * vy + vz * vz);

	  mesh[i][0] _MASS = rho;
	  mesh[i][0] _MOMX = rho * vx;
	  mesh[i][0] _MOMY = rho * vy;
	  mesh[i][0] _MOMZ = rho * vz;
	  mesh[i][0] _ENER = pressure * gammam1i + bsqr + ke;
	  mesh[i][0] _B_X = bx;
	  mesh[i][0] _B_Y = by;
	  mesh[i][0] _B_Z = bz;
	  mesh[i][0].temperature = myid;

	  mesh[i][1] = mesh[i][0];

	}
    }

  if (myWest == MPI_PROC_NULL)
    {
      // Lower xboundary is reflecting axisymmetric
      for (int i = 0; i < ny - 0; ++i)
	{
	  recvmesh[1][i] = mesh[2][i];
	  recvmesh[0][i] = mesh[3][i];
#ifdef CYLINDRICAL
	  recvmesh[0][i] _MOMX = (-mesh[0][i] _MOMX);
	  recvmesh[1][i] _MOMX = (-mesh[1][i] _MOMX);
	  recvmesh[0][i] _MOMZ = (-mesh[0][i] _MOMZ);
	  recvmesh[1][i] _MOMZ = (-mesh[1][i] _MOMZ);
	  recvmesh[0][i] _B_X = (-mesh[0][i] _B_X);
	  recvmesh[1][i] _B_X = (-mesh[1][i] _B_X);
	  recvmesh[0][i] _B_Z = (-mesh[0][i] _B_Z);
	  recvmesh[1][i] _B_Z = (-mesh[1][i] _B_Z);
#endif
	}
    }


  if (myEast == MPI_PROC_NULL)
    {
      // Upper xboundary is freeflow
      for (int j = 0; j < ny; ++j)
	{
	  recvmesh[nx - 2][j] = mesh[nx - 3][j];
	  recvmesh[nx - 1][j] = mesh[nx - 2][j];
	}
    }


  if (myNorth == MPI_PROC_NULL)
    {
      // Upper xboundary is freeflow
      for (int i = 0; i < nx; ++i)
	{
	  recvmesh[i][ny - 2] = mesh[i][ny - 3];
	  recvmesh[i][ny - 1] = mesh[i][ny - 2];
	}
    }

#endif

	std::cout << __FUNCTION__ << " By 202,2 = " << mesh[202][2] _B_Y << std::endl; 
  return 0;
}
