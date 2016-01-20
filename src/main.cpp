/* $Id: main.cpp,v 1.6 2006-11-16 13:48:07 gmurphy Exp $  */
/* 
 * PHYSICS:
 * Ideal MHD - Cartesian
 * Resistivity
 * Axisymmetry
 * Gravity
 * Staggered Mesh
 * Limiter
 *
 *
 * BUGS:
 * Axisymmetric Boundary conditions /reflecting 
 *
 * 
 * DESCRIPTION:
 * 
 * MPI program using of derived datatypes to send boundary
 * conditions. With this system can send a block of boundary conditions to
 * another processor. The datatype used is a
 * Template Numerical Toolkit 2D Array of structs
 *
 * FUTURE WORK:
 *
 * BOUNDARY:
 * Boundary cells send to/receive from  MPI_PROC_NULL
 *
 * NOTES: 
 * Am using MPI_UB to make sure the data is the correct size.
 *
 * */
#include <mpi.h>
#include <tnt.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#define NX 160
#define NY 160

#include "maxspeed.h"
#include "out.h"
#include "initialise_maes.h"
#include "initialise_blast.h"
#include "initialise_uniform.h"
#include "orszagtang.h"
#include "zanni.h"
#include "zannisimple.h"
#include "initialise_jet.h"
#include "parboundary.h"
#include "parflux.h"
#include "parupdate.h"
#include "problem.h"
#include "physics.h"
#include "emf_exchange.h"

#define ROOT 0
int
main (int argc, char *argv[])
{
  int myid = 99, numprocs;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int namelen;
  // MPI kickoff
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name (processor_name, &namelen);

  std::ifstream init;

  double cfl = 0.40;
  int nx = NX;
  int ny = NY;
  int ne = NE;
  int printtime = 10;
  int maxstep = 5000;


  if (myid == ROOT)
    {
      init.open ("setup.dat");
      init >> printtime;
      init.ignore (256, '\n');
      init >> cfl;
      init.ignore (256, '\n');
      init >> nx;
      init.ignore (256, '\n');
      init >> ny;
      init.ignore (256, '\n');
      init >> maxstep;
      init.ignore (256, '\n');
      init.close ();
    }

  MPI_Bcast (&printtime, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&cfl, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&maxstep, 1, MPI_INT, 0, MPI_COMM_WORLD);

  double rgi = 0;
  double gammag = PhysConsts::gamma;
  double gammam1 = gammag - 1.;
  double gammam1i = 1. / gammam1;
  double delta_x = 0.001;
  double delta_y = delta_x;
  TNT::Array2D < unk > mesh (nx, ny);
  TNT::Array2D < unk > halfmesh (nx, ny);
  TNT::Array2D < unk > newmesh (nx, ny);

  //                                                            (nx-1 , ny)
  TNT::Array2D < flux > fx (nx , ny);
  //                                                            (nx , ny-1)
  TNT::Array2D < flux > fy (nx, ny);
  TNT::Array2D < flux > SecondOrdCorr (nx, ny);

  TNT::Array2D < double >faceBx (nx + 1, ny + 1);
  TNT::Array2D < double >faceBy (nx + 1, ny + 1);
  TNT::Array2D < double >faceBxt (nx + 1, ny + 1);
  TNT::Array2D < double >faceByt (nx + 1, ny + 1);
  TNT::Array2D < double >emfx (nx , ny );
  TNT::Array2D < double >emfy (nx , ny );
  TNT::Array1D < double >SoundSpeedMidPlane (nx);

  // Current arrays
  TNT::Array2D < double >xjx (nx, ny);
  TNT::Array2D < double >xjy (nx, ny);
  TNT::Array2D < double >xjz (nx, ny);
  TNT::Array2D < double >yjx (nx, ny);
  TNT::Array2D < double >yjy (nx, ny);
  TNT::Array2D < double >yjz (nx, ny);
  emfx = 0;
  emfy = 0;

  int rank_above = 0, rank_below = 0, tag = 0;
  MPI_Datatype MpiColumnType;
  MPI_Datatype MpiRowType;
  MPI_Datatype MPI_unk_type;
  int ii = 0, jj = 0;
  MPI_Status status;
  unk test;
  unk maxvar;
  unk minvar;
  //MPI_Aint disps[3] = { 0, (int) test.array - (int) &test.temperature, sizeof (test) };	/* guessing... */
  MPI_Aint disps[3] ;
  MPI_Aint addr[3] ;
  //MPI_Get_Address (test.array, &addr[0]  );
  disps[0]= 0;
  disps[1]=sizeof (double);
  //disps[1]= test.array -  &test.temperature;
  disps[2]=sizeof (test) ;	/* guessing... */
  int blks[3] = { 1, ne, 1 };
  MPI_Datatype types[3] = { MPI_INT, MPI_DOUBLE, MPI_UB };
  int testint = 0;
  int tstep = 0;
  std::ofstream outfile;
  int rc = 0;
  double maxspeed_ = 0;
  double tmpbuf = 0;
  double pressure = 0;
  double rr, px, py, pz, et, ri, vx, vy, vz, ke, p, ptot;
  double bx, by, bz, b2, vdotb;
  MPI_Comm old_comm, Cart_comm;
  int ndims, reorder, ierr;
  int dim_size[2], periods[2];

  int arr_size[2];
  double time = 0.0;
  const double eps = 0.1;
  arr_size[0] = nx;
  arr_size[1] = ny;
  TNT::Array2D < double >smallarr (nx, ny);





  /*
     printf ("Process %d on %s\n", myid, processor_name);
     MPI_Barrier (MPI_COMM_WORLD);
   */


  MPI_Type_struct (3, blks, disps, types, &MPI_unk_type);
  MPI_Type_commit (&MPI_unk_type);


  old_comm = MPI_COMM_WORLD;
  ndims = 2;
  dim_size[0] = 0;
  dim_size[1] = 0;
  periods[0] = 0;
  periods[1] = 0;
#ifdef   PERIODIC
  periods[0] = 1;
  periods[1] = 1;
#endif
  reorder = 1;



  MPI_Dims_create (numprocs, ndims, dim_size);
  if (dim_size[0] > dim_size[1])
    {
      int a = dim_size[0];
      dim_size[0] = dim_size[1];
      dim_size[1] = a;
    }

  // Could substitute own dim creation function

  MPI_Cart_create (old_comm, ndims, dim_size, periods, reorder, &Cart_comm);


  int coords[] = { 0, 0 };

  MPI_Cart_coords (Cart_comm, myid, ndims, coords);


  int myaddress[] = { 0, 0 };
  myaddress[0] = (nx - 2 * NGC) * coords[0];
  myaddress[1] = (ny - 2 * NGC) * coords[1];



  if (myid == ROOT)
    {
      init.open ("setup.dat");
      init >> printtime;
      init.ignore (256, '\n');
      init.close ();
      system ("rm -rf output/*.log");


      std::cout << "\t\t\t2.5D MAES Solver Code" << std::endl;
      std::cout << "\t\t\tVersion 2.0" << std::endl;
      std::
	cout << "\t\t\tCompiled on " << __DATE__ << " at " << __TIME__ <<
	std::endl;

      std::cout << "\t\t\tNo of steps 	=  " << maxstep << std::endl;
      std::cout << "\t\t\tnx  		 	=  " << nx << std::endl;
      std::cout << "\t\t\tny    			=  " << ny << std::
	endl;
      std::cout << "\t\t\tne    			=  " << ne << std::
	endl;
      std::cout << "\t\t\tCFL   			=  " << cfl << std::
	endl;
      std::
	cout << "\t\t\tGamma 			=  " << gammag << std::endl;
      std::cout << "\t\t\tPrint time	=  " << printtime << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
    }


  printf ("Process %d on %s\n", myid, processor_name);
  MPI_Barrier (MPI_COMM_WORLD);

  float myx = 0;
  float myy = 0;

// Initialisation
  float test_float;
  MPI_Bcast (&test_float, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  unk JetVal;
  if (myid == ROOT)
    {
      JetVal _MASS = 999999;
    }
  MPI_Bcast (&JetVal, 1, MPI_unk_type, 0, MPI_COMM_WORLD);
  std::cout << JetVal _MASS << std::endl;

  double xsize = 40;
  double ysize = 80;
#ifdef MAES
  //initialise_maes (mesh, faceBx, faceBy, Cart_comm, ndims, dim_size, &xsize, &ysize);
  //zanni (mesh, faceBx, faceBy, Cart_comm, ndims, dim_size, &xsize, &ysize);
  zannisimple (mesh, faceBx, faceBy, Cart_comm, ndims, dim_size, &xsize, &ysize);
#endif

#ifdef BLAST
  initialise_blast (mesh,
		    faceBx,
		    faceBy, Cart_comm, ndims, dim_size, &xsize, &ysize);
#endif
#ifdef ORSZAGTANG
  orszagtang (mesh,
		      faceBx,
		      faceBy, Cart_comm, ndims, dim_size, &xsize, &ysize);
#endif
#ifdef JET
  initialise_jet (mesh,
		  faceBx, faceBy, Cart_comm, ndims, dim_size, &xsize, &ysize);
#endif

  newmesh = mesh.copy ();
  delta_x = xsize / (dim_size[0] * (nx - 2 * NGC));
  delta_y = ysize / (dim_size[1] * (ny - 2 * NGC));

  newmesh = mesh.copy ();
  delta_x = xsize / (dim_size[0] * (nx - 2 * NGC));
  delta_y = ysize / (dim_size[1] * (ny - 2 * NGC));

#ifdef DEBUG_OUT
  std::cout << "Proc" << myid << " Good to here: " << __LINE__ << std::endl;
  MPI_Barrier (MPI_COMM_WORLD);
#endif



#ifdef DEBUG_OUT
  std::cout << "Proc" << myid << " Good to here: " << __LINE__ << std::endl;
  MPI_Barrier (MPI_COMM_WORLD);
#endif
  // out (numprocs, myid, Cart_comm, ndims, dim_size, mesh, arr_size);


  MPI_Type_vector (nx,		/* # column elements */
		   1,		/* 1 column only */
		   ny,		/* skip ny elements */
		   MPI_unk_type,	/* elements are double */
		   &MpiColumnType);	/* MPI derived datatype */

  MPI_Type_vector (1,		/* # row elements */
		   ny,		/* 1 row only */
		   0,		/* skip ny elements */
		   MPI_unk_type,	/* elements are double */
		   &MpiRowType);	/* MPI derived datatype */

  MPI_Type_commit (&MpiColumnType);
  MPI_Type_commit (&MpiRowType);

  /*
   *    Locate my neighbours, by using my coords
   *
   *    My North Neighbour = coords[1]+1
   *    My South Neighbour = coords[1]-1
   *    My East Neighbour = coords[0]+1
   *    My West Neighbour = coords[0]-1
   * */
  int myNorth = 0;
  int mySouth = 0;
  int myEast = 0;
  int myWest = 0;

  int tcoords[] = { 0, 0 };
  tcoords[0] = coords[0];
  tcoords[1] = coords[1];
#ifndef PERIODIC
  tcoords[0] = coords[0];
  tcoords[1] = coords[1] + 1;
  if (tcoords[1] >= dim_size[1])
    {
      myNorth = MPI_PROC_NULL;
    }
  else
    {
      MPI_Cart_rank (Cart_comm, tcoords, &myNorth);
    }
  tcoords[0] = coords[0];
  tcoords[1] = coords[1] - 1;
  if (tcoords[1] < 0)
    {
      mySouth = MPI_PROC_NULL;
    }
  else
    {
      MPI_Cart_rank (Cart_comm, tcoords, &mySouth);
    }
  tcoords[0] = coords[0] + 1;
  tcoords[1] = coords[1];
  if (tcoords[0] >= dim_size[0])
    {
      myEast = MPI_PROC_NULL;
    }
  else
    {
      MPI_Cart_rank (Cart_comm, tcoords, &myEast);
    }
  tcoords[0] = coords[0] - 1;
  tcoords[1] = coords[1];
  if (tcoords[0] < 0)
    {
      myWest = MPI_PROC_NULL;
    }
  else
    {
      MPI_Cart_rank (Cart_comm, tcoords, &myWest);
    }
#else
  tcoords[0] = coords[0];
  tcoords[1] = coords[1] + 1;
  MPI_Cart_rank (Cart_comm, tcoords, &myNorth);
  tcoords[0] = coords[0];
  tcoords[1] = coords[1] - 1;
  MPI_Cart_rank (Cart_comm, tcoords, &mySouth);
  tcoords[0] = coords[0] + 1;
  tcoords[1] = coords[1];
  MPI_Cart_rank (Cart_comm, tcoords, &myEast);
  tcoords[0] = coords[0] - 1;
  tcoords[1] = coords[1];
  MPI_Cart_rank (Cart_comm, tcoords, &myWest);
#endif


  rank_above = myid + 1;
  if (rank_above > numprocs - 1)
    {
      rank_above = 0;
      rank_above = MPI_PROC_NULL;
    }
  rank_below = myid - 1;
  if (rank_below < 0)
    {
      rank_below = numprocs - 1;
      rank_below = MPI_PROC_NULL;
    }

  std::cout
    << " MyId " << myid
    << " MyEast " << myEast
    << " MyNorth " << myNorth
    << " MySouth " << mySouth << " MyWest " << myWest << std::endl;




#define DEBUG_OUT3
#ifdef DEBUG_OUT3
  std::cout << "Proc" << myid << " Good to here: " << __LINE__ << std::endl;
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  outhdf5 (numprocs, myid, Cart_comm, ndims, dim_size, mesh, arr_size, 7777);

  parboundary (mesh, mesh,  faceBx, faceBy, myNorth, mySouth, myEast, myWest, myid,
	       Cart_comm);

  faceBxt = faceBx.copy ();
  faceByt = faceBy.copy ();
  outhdf5 (numprocs, myid, Cart_comm, ndims, dim_size, mesh, arr_size, tstep);

  TNT::Array2D < double >divb (nx, ny);
  double idx = 1.0 / delta_x;
  double idy = 1.0 / delta_y;




  // Loop in time
  for (tstep = 1; tstep < maxstep; tstep++)
    {



#define CHECK_DIVB
#ifdef CHECK_DIVB
      double maxdivb = (-9e99);
      double mindivb = 9e99;
      for (int i = 2; i < nx - 2; i++)
	for (int j = 2; j < ny - 2; j++)
	  {

	    divb[i][j] =
	      idx * (faceBx[i + 1][j] - faceBx[i][j]) +
	      idy * (faceBy[i][j + 1] - faceBy[i][j]);

#ifdef CYLINDRICAL
	    double myx = (double) myaddress[0] + (i - NGC);
	    double dri = 2.0 / (2 * myx + 1);
	    divb[i][j] = idx * dri * ((myx + 1.0) * faceBx[i + 1][j] - (myx ) * faceBx[i][j]) + idy * (faceBy[i][j + 1] - faceBy[i][j]);
	    //divb[i][j] = idx * dri * ((myx + 0.5) * faceBx[i + 1][j] - (myx -0.5) * faceBx[i][j]) + idy * (faceBy[i][j + 1] - faceBy[i][j]);
#endif


	    maxdivb = std::max (maxdivb, divb[i][j]);
	    mindivb = std::min (mindivb, divb[i][j]);


	  }

      std::cout << std::setiosflags (std::ios::scientific);
      std::cout << __FUNCTION__ << std::endl;
      std::cout << " MaxDivB = " << maxdivb << std::endl;
      std::cout << " MinDivB = " << mindivb << std::endl;


#endif



#ifdef STAGGER_MESH
      // Average face-centred B fields to cell centres
      for (ii = 2; ii < nx - 2; ii++)
	for (jj = 2; jj < ny - 2; jj++)
	  {
	    mesh[ii][jj] _B_X = 0.5 * (faceBx[ii + 1][jj] + faceBx[ii][jj]);
#ifdef CYLINDRICAL
	    //double myx = myaddress[0] + (ii - NGC);
	    //rgi = (myx * myx + myx + 1.0 / 3.0) * delta_x / (myx + 0.5);
	   //gg rgi = (myx * myx - myx + 1.0 / 3.0) * delta_x / (myx - 0.5);
	    //mesh[ii][jj] _B_X = faceBx[ii + 1][jj] - idx * ((ii + 1) * delta_x - rgi) * (faceBx[ii + 1][jj] - faceBx[ii][jj]);
	    mesh[ii][jj] _B_X = 0.5 * (faceBx[ii + 1][jj] + faceBx[ii][jj]);
#endif
	    mesh[ii][jj] _B_Y = 0.5 * (faceBy[ii][jj + 1] + faceBy[ii][jj]);
	  }
      parboundary (mesh, mesh, faceBx, faceBy, myNorth, mySouth, myEast, myWest,
		   myid, Cart_comm);
#endif

      if (myid == ROOT)
	{
	  std::cout << "Timestep = " << tstep;
	}

      for (int qq = 0; qq < 8; ++qq)
	{
	  maxvar.array[qq] = (-9e99);
	  minvar.array[qq] = 9e99;
	}
      // time Advection

      maxspeed_ = 0.0;
      // Compute the maximum wave speed
      rc = maxspeed (mesh, mesh[ii][jj].array, &maxspeed_);
      if (rc == 1)
	{
	  std::cout << __FUNCTION__ << ": "
	    << " proc: " << myid
	    << " (" << ii << " ," << jj << ")" << std::endl;
	  MPI_Abort (MPI_COMM_WORLD, 1);
	  exit (0);
	}


      // Communicate the max speed to each processor
      MPI_Barrier (MPI_COMM_WORLD);
      MPI_Allreduce (&maxspeed_,
		     &tmpbuf, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      maxspeed_ = tmpbuf;
      if (myid == ROOT)
	{
	  //std::cout << "Max Speed = " << maxspeed_ << std::endl;
	  if (std::isinf (maxspeed_))
	    {
	      std::
		cout << "Max is  " << " maxspeed = " << maxspeed_ <<
		" , dx = " << delta_x << std::endl;
	      MPI_Abort (MPI_COMM_WORLD, 7777);
	      exit (0);
	    }
	}

      double del;

      del = cfl / maxspeed_;
      //del=del*del;
      double delta_t = del * delta_x;
      double dtodx = del;
      time = time + delta_t;

      if (std::isnan (del))
	{
	  std::
	    cout << "Del is nan " << " maxspeed = " << maxspeed_ << " , dx = "
	    << delta_x << std::endl;
	  MPI_Abort (MPI_COMM_WORLD, 7777);
	  exit (0);
	}
      if (std::isinf (del))
	{
	  std::cout << "Del is inf "
	    << " maxspeed = " << maxspeed_ << " , dx = " << delta_x << std::
	    endl;
	  MPI_Abort (MPI_COMM_WORLD, 7777);
	  exit (0);
	}
      if (myid == ROOT)
	{
	  std::cout << "\tTime= " << std::setiosflags (std::ios::
						       scientific) << time;
	  std::cout << "\tMaxspeed= " << std::setiosflags (std::ios::
							   scientific) <<
	    maxspeed_;
	  std::cout 
	  << "\tCFL= " 
	  << std::setiosflags (std::ios:: scientific)
	  << del;
	  std::cout << std::endl;


	}



      // Compute x and y fluxes
      parflux (mesh, fx, fy, xjx, xjy, xjz, yjx, yjy, yjz, emfx, emfy, faceBx, faceBy, SoundSpeedMidPlane, del, delta_x, delta_y, myaddress, tstep, Cart_comm);
      //parsecond (mesh, fx, fy, SecondOrdCorr, xjx, xjy, xjz, yjx, yjy, yjz, emfx, emfy, faceBx, faceBy, SoundSpeedMidPlane, del, delta_x, delta_y, myaddress, Cart_comm);
      // communicate emf (from the fluxes) between processors
      emf_exchange (emfx, emfy, myNorth, mySouth, myEast, myWest, myid);

//#define SECOND_ORDER_STEP
#ifdef SECOND_ORDER_STEP
      // Do time update onto newmesh
      parupdate (newmesh, mesh,  mesh, fx, fy, SecondOrdCorr, faceBxt, faceByt, xjx, xjy, xjz, yjx, yjy, yjz, emfx, emfy, SoundSpeedMidPlane, 0.5*delta_t, delta_x, delta_y, Cart_comm, 0);

//    get rid of me !!!!!   parboundary (newmesh, newmesh, faceBxt, faceByt, myNorth, mySouth, myEast, myWest, myid, Cart_comm);

      // Copy newmesh onto oldmesh
      halfmesh = newmesh.copy ();

#ifdef STAGGER_MESH
         // Average face-centred B fields to cell centres
         for (ii = 2; ii < nx - 2; ii++)
         for (jj = 2; jj < ny - 2; jj++)
         {
         halfmesh[ii][jj] _B_X = 0.5 * (faceBxt[ii + 1][jj] + faceBxt[ii][jj]);
#ifdef CYLINDRICAL
	    //double myx = myaddress[0] + (ii - NGC);
	    //rgi = (myx * myx + myx + 1.0 / 3.0) * delta_x / (myx + 0.5);
	    //gg rgi = (myx * myx - myx + 1.0 / 3.0) * delta_x / (myx - 0.5);
	    //halfmesh[ii][jj] _B_X = faceBxt[ii + 1][jj] - idx * ((ii + 1) * delta_x - rgi) * (faceBxt[ii + 1][jj] - faceBxt[ii][jj]);
         halfmesh[ii][jj] _B_X = 0.5 * (faceBxt[ii + 1][jj] + faceBxt[ii][jj]);
#endif

         halfmesh[ii][jj] _B_Y = 0.5 * (faceByt[ii][jj + 1] + faceByt[ii][jj]);
         }
      // get rid of me!!!  parboundary (halfmesh, mesh,  faceBxt, faceByt, myNorth, mySouth, myEast, myWest, myid, Cart_comm);
#endif
      parboundary (halfmesh, halfmesh,  faceBxt, faceByt, myNorth, mySouth, myEast, myWest, myid, Cart_comm);

      MPI_Barrier (MPI_COMM_WORLD);
      if (myid == ROOT) { std::cout << "Halfstep" << std::endl; }
      // Compute x and y fluxes
      parsecond (halfmesh, fx, fy, SecondOrdCorr, xjx, xjy, xjz, yjx, yjy, yjz, emfx, emfy, faceBxt, faceByt, SoundSpeedMidPlane, del, delta_x, delta_y, myaddress, Cart_comm);
      // communicate emfs between processors
      emf_exchange (emfx, emfy, myNorth, mySouth, myEast, myWest, myid);
      // Do time update onto newmesh
      parupdate (newmesh, halfmesh, mesh, fx, fy, SecondOrdCorr, faceBx, faceBy, xjx,
		 xjy, xjz, yjx, yjy, yjz, emfx, emfy, SoundSpeedMidPlane,
		 delta_t, delta_x, delta_y, Cart_comm, 1);
#else
      if (myid == ROOT) { std::cout << "First order" << std::endl; }
      parupdate (newmesh, mesh, mesh, fx, fy, SecondOrdCorr, faceBx, faceBy, xjx, xjy, xjz, yjx, yjy, yjz, emfx, emfy, SoundSpeedMidPlane, delta_t, delta_x, delta_y, Cart_comm, 0);
#endif

#ifdef DEBUGG
   std::cout << __FUNCTION__ << " By 201,2 = " << mesh[202][3] _B_Y << std::endl;    
   std::cout << __FUNCTION__ << " By 202,2 = " << mesh[202][2] _B_Y << std::endl;    
#endif


      // Copy newmesh onto oldmesh
      mesh = newmesh.copy ();




      // Do boundaries
      parboundary (mesh, mesh,  faceBx, faceBy, myNorth, mySouth, myEast, myWest,
		   myid, Cart_comm);
   std::cout << __FUNCTION__ << " By 202,2 = " << mesh[202][2] _B_Y << std::endl;    
      faceBxt = faceBx.copy ();
      faceByt = faceBy.copy ();

      for (ii = 2; ii < nx - 2; ii++)
	for (jj = 2; jj < ny - 2; jj++)
	  {


	    for (int k = 0; k < ne; k++)
	      {

		maxvar.array[k] =
		  (maxvar.array[k] >
		   mesh[ii][jj].array[k] ? maxvar.
		   array[k] : mesh[ii][jj].array[k]);
		minvar.array[k] =
		  (minvar.array[k] <
		   mesh[ii][jj].array[k] ? minvar.
		   array[k] : mesh[ii][jj].array[k]);
	      }
#ifdef CHECK_EN
	    if (minvar _ENER < 0)
	      {
		std::cout
		  << "Negative energy "
		  << minvar _ENER
		  << " Location: (" << ii << "," << jj << ")" << std::endl;
		MPI_Abort (MPI_COMM_WORLD, 5605050);
		exit (0);
	      }
#endif

	  }


      unk maxv;
      unk minv;

      for (int qq = 0; qq < ne; qq++)
	{
	  MPI_Reduce (&maxvar.array[qq], &maxv.array[qq],
		      1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	  MPI_Reduce (&minvar.array[qq], &minv.array[qq],
		      1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	}

      if (myid == ROOT)
	{
	  for (int k = 0; k < ne; k++)
	    {
//	      std::cout << std::setiosflags (std::ios::fixed) << std::setprecision(10);
	      std::cout << k +
		1 << " " << maxv.array[k] << " " << minv.
		array[k] << std::endl;
	    }
	}

      // Output
      if (tstep % printtime == 0 || tstep < 10)
	{



	  TNT::Array2D < unk > outmesh (nx, ny);
	  outmesh = mesh.copy ();



#ifdef STAGGER_MESH
	  // Average face-centred B fields to cell centres
	  for (ii = 2; ii < nx - 2; ii++)
	    for (jj = 2; jj < ny - 2; jj++)
	      {
		outmesh[ii][jj] _B_X = 0.5 * (faceBx[ii + 1][jj] + faceBx[ii][jj]);
#ifdef CYLINDRICAL
		//double myx = myaddress[0] + (ii - NGC);
		//rgi = (myx * myx + myx + 1.0 / 3.0) * delta_x / (myx + 0.5);
		//gg rgi = (myx * myx - myx + 1.0 / 3.0) * delta_x / (myx - 0.5);
		//outmesh[ii][jj] _B_X = faceBx[ii + 1][jj] - idx * ((ii + 1) * delta_x - rgi) * (faceBx[ii + 1][jj] - faceBx[ii][jj]);
		outmesh[ii][jj] _B_X = 0.5 * (faceBx[ii + 1][jj] + faceBx[ii][jj]);
#endif
		outmesh[ii][jj] _B_Y =
		  0.5 * (faceBy[ii][jj + 1] + faceBy[ii][jj]);
	      }
	  parboundary (outmesh, outmesh, faceBx, faceBy, myNorth, mySouth, myEast, myWest, myid, Cart_comm);
#endif

	if (myid == ROOT)
	{
		std::ofstream opf;
		opf.open("xxx.txt");
	  for (ii = 0; ii < nx - 1; ii++)
		  {
			 opf << (mesh [ii][20] _MASS) << std::endl; 
			  }
	}
	  outhdf5 (numprocs, myid, Cart_comm, ndims, dim_size, outmesh, arr_size,
	       tstep);
	  MPI_Barrier (MPI_COMM_WORLD);
	}
      // Gather all of small arrays into a large array and print out. 

    }



  // Final Output
  outhdf5 (numprocs, myid, Cart_comm, ndims, dim_size, mesh, arr_size, tstep);


  MPI_Type_free (&MPI_unk_type);
  MPI_Type_free (&MpiColumnType);
  MPI_Finalize ();
  std::cout 
  << " Done "
  <<std::endl;
  return 0;
}
