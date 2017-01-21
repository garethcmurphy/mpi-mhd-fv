/*
 *
 */


#include <mpi.h>
#include <fstream>
#include <sstream>
#include "tnt.h"
#include "out.h"

//using namespace std;
//using namespace TNT;

#define NE 9
#define SCALE 1




int
out (int numprocs, int myid, MPI_Comm new_comm, int ndims, int *dim_size,
     TNT::Array2D < unk > small_mesh, int *arr_size, int step)
{
  const int ne = NE;
  const int small_x = arr_size[0];
  const int small_y = arr_size[1];
  TNT::Array2D < unk > spare_mesh (200, 200);
  TNT::Array2D < double >small_arr (small_x, small_y);
  int big_x = 9;
  int big_y = 9;
  big_x = (small_x - 2 * NGC) * (dim_size[0]);
  big_y = (small_y - 2 * NGC) * (dim_size[1]);
  MPI_Datatype big_type;
  MPI_Datatype small_row;
  MPI_Datatype small_block;
  MPI_Status status;
  int tag = 0;

  double big_arr[big_x][big_y];
  TNT::Array2D < unk > big_mesh (big_x, big_y);

  int *disp;
  int *fiend;
  int *size;
  int *order;
  int coords[2];
  std::ofstream outFile;




  MPI_Datatype mycolumn;

  unk test;
  //MPI_Aint displacements[3] = { 0, (int) test.array - (int) &test.temperature, sizeof (test) };	/* guessing... */
    MPI_Aint displacements[2];
    displacements[0] = sizeof(double);
    displacements[1] = sizeof(test);
    int blks[2] = {1, ne};
    MPI_Datatype types_str[2] = {MPI_INT, MPI_DOUBLE};


    MPI_Type_create_struct (3, blks, displacements, types_str, &mycolumn);
  MPI_Type_commit (&mycolumn);



  fiend = new int[numprocs];
  order = new int[numprocs];
  size = new int[numprocs];
  disp = new int[numprocs];


  // Initialise big_array
  for (int j = 0; j < big_y; j++)
    for (int i = 0; i < big_x; i++)
      {
	big_arr[i][j] = 99;
	big_mesh[i][j].temperature = 99;
      }

  // Initialise small_arr
  for (int j = 0; j < small_y; j++)
    for (int i = 0; i < small_x; i++)
      {
	small_arr[i][j] = myid * SCALE;
	//    small_mesh[i][j].temperature = myid*SCALE;
      }



//  std::cout << std::endl;
//  std::cout << std::endl;

  for (int j = 0; j < numprocs; j++)
    {
      size[j] = 1;
      order[j] = j;
    }

  // Grab the coords of the processors and send them all to the root processor
  MPI_Cart_coords (new_comm, myid, ndims, coords);
  int tmp = coords[0] + coords[1] * big_x;
//  std::cout << coords[0] << " " << coords[1] << " " << tmp << std::endl;
  MPI_Gatherv (&tmp, 1, MPI_INT, disp, size, order, MPI_INT, 0, new_comm);

  MPI_Gatherv (&myid, 1, MPI_INT, fiend, size, order, MPI_INT, 0, new_comm);

  // Define a row vector

    MPI_Aint sizeofdouble = 0;
    MPI_Aint sizeofunk = 0;
    MPI_Aint lb;
    MPI_Type_get_extent(MPI_DOUBLE, &lb, &sizeofdouble);
    MPI_Type_get_extent(mycolumn, &lb, &sizeofunk);

  // Define a datatype with a stride 
//  MPI_Type_vector (small_x, 1, 2, MPI_DOUBLE, &small_row);
//  MPI_Type_hvector (small_y, 1, big_x*sizeofdouble, MPI_DOUBLE, &small_block);

    MPI_Type_create_hvector((small_x - 2 * NGC), (small_y - 2 * NGC),
		    big_y * sizeofunk, mycolumn, &small_row);
  MPI_Type_commit (&small_row);

    MPI_Aint disps[1] = {(small_x - 2 * NGC) * sizeof(unk)};    /* guessing... */
    int blocklengths[1] = {1};
    MPI_Datatype types[1] = {small_row};

    MPI_Type_create_struct(1, blocklengths, disps, types, &big_type);
  MPI_Type_commit (&big_type);

//  MPI_Gatherv (&small_arr[0][0], small_x * small_y, MPI_DOUBLE,
//             big_arr, size, disp, big_type, 0, new_comm);

  MPI_Type_vector ((small_x - 2 * NGC), (small_y - 2 * NGC), small_x,
		   mycolumn, &small_block);
  MPI_Type_commit (&small_block);

#undef DEBUG_OUT
#ifdef DEBUG_OUT
  std::
    cout << "Proc" << myid << " Good to before exchange: " << __FUNCTION__ <<
    " " << __LINE__ << std::endl;
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  if (myid == 0)
    {
      //std::cout << "Large Mesh " << big_x << "," << big_y << std::endl;
      for (int i = 1; i < numprocs; i++)
	{
	  MPI_Cart_coords (new_comm, i, ndims, coords);
#ifndef DEBUGG

	  MPI_Recv (&big_mesh[(small_x - 2 * NGC) * coords[0]]
		    [(small_y - 2 * NGC) * coords[1]], 1, big_type, i, i,
		    MPI_COMM_WORLD, &status);
#else
	  std::cout
	    << "Receiving small mesh"
	    << "(" << small_x << "," << small_y << ")"
	    << "into coords ("
	    << (small_x - 2 * NGC) * coords[0]
	    << " , " << (small_y - 2 * NGC) * coords[1] << ") " << std::endl;
	  MPI_Recv (&spare_mesh[0][0], 1, big_type, i, i, MPI_COMM_WORLD,
		    &status);
#endif

	}
      for (int j = 0; j < small_y - 2 * NGC; j++)
	for (int i = 0; i < small_x - 2 * NGC; i++)
	  {
	    //    big_arr[i][j] = small_arr[i][j];
	    big_mesh[i][j] = small_mesh[i + NGC][j + NGC];
	  }
    }
  else
    {
      MPI_Send (&small_mesh[NGC][NGC], 1, small_block, 0, myid,
		MPI_COMM_WORLD);
    }

  if (myid == 0)
    {
      /*
         std::cout << "Root proc" << std::endl;
         std::cout << dim_size[0] << " " << dim_size[1] << std::endl;
         for (int j = 0; j < numprocs; j++)
         std::cout << "\t " << fiend[j];
         std::cout << std::endl;
         std::cout << std::endl;
         for (int j = 0; j < numprocs; j++)
         std::cout << "\t " << disp[j];
         std::cout << std::endl;
         std::cout << std::endl;
       */


      std::string filename = "out_";
      std::string blockname = "out_";
      std::stringstream stream_filename;
      std::stringstream stream_tag;
      std::string outputdir = "./";
      std::string str_file_tag;
      stream_tag.clear ();
      stream_tag.width (4);
      stream_tag.fill ('0');
      stream_tag << step;
      stream_tag >> str_file_tag;
      stream_filename.clear ();
      stream_filename << outputdir << blockname << str_file_tag << ".cut";
      stream_filename >> blockname;


#ifdef PGM_FILE
      outFile.open (blockname.c_str ());
      outFile << "P2" << std::endl;;
      outFile << "# outFile.pgm " << std::endl;;
      outFile << big_x << " ";
      outFile << big_y << std::endl;
      outFile << numprocs << std::endl;

      for (int j = 0; j < big_y; j++)
	{
	  for (int i = 0; i < big_x; i++)
	    {
	      outFile << " " << big_mesh[i][j].temperature;
	    }
	  outFile << std::endl;
	}
      outFile << std::endl;

      outFile.close ();
#endif

      stream_filename.clear ();
      stream_filename << outputdir << filename << str_file_tag << ".dat";
      stream_filename >> filename;

      if (myid == 0)
	{
	  std::cout << "Proc. " << myid << " ecrit a fichier ";
	  std::cout << filename << std::endl;
	}

  TNT::Array2D < double  >  az (big_x, big_y);
      for (int j = 1; j < big_y; j++)
	{
	  for (int i = 1; i < big_x; i++)
	    {
			 
			 double dbydx=0.;
			 double dbxdy=0.;
			 dbydx=  big_mesh [i][j] _B_Y - big_mesh [i-1][j] _B_Y;
			 dbxdy=  big_mesh [i][j] _B_X - big_mesh [i][j-1] _B_X;
			 az [i][j]=  dbydx  - dbxdy; 

		 }
	}

      //outFile.open ("density.dat");
      outFile.open (filename.c_str ());
      outFile << "# " << big_x << std::endl;
      outFile << "# " << big_y << std::endl;
      for (int j = 0; j < big_y; j++)
	{
	  for (int i = 0; i < big_x; i++)
	    {
//#define CHECKSYM

#ifdef CHECKSYM
			double a =big_mesh[i][j] _MASS;
			double b =  big_mesh[big_x-i-1][big_y-j-1] _MASS;
			if ( a - b > 1e-10)
			{
			std::cout 
         << i << " " << j
         << " " << big_x-i-1 << " " << big_y-j-1 
			<< " " << a <<  " " << b 
			<< " " << a - b 
			<< std::endl;
			}
#endif

	      outFile << " " << i;
	      outFile << " " << j;
	      outFile << " " << log10 (big_mesh[i][j] _MASS);
	      for (int qq = 0; qq < NE; qq++)
		{
		  outFile << " " << big_mesh[i][j].array[qq];
		}

	      double gammam1 = PhysConsts::gamma - 1;
	      double rr = big_mesh[i][j] _MASS;
	      double px = big_mesh[i][j] _MOMX;
	      double py = big_mesh[i][j] _MOMY;
	      double pz = big_mesh[i][j] _MOMZ;
	      double et = big_mesh[i][j] _ENER;
	      double bx = big_mesh[i][j] _B_X;
	      double by = big_mesh[i][j] _B_Y;
	      double bz = big_mesh[i][j] _B_Z;
	      double ri = 1.0 / rr;
	      double vx = px * ri;
	      double vy = py * ri;
	      double vz = pz * ri;
	      double ke = 0.5 * rr * (vx * vx + vy * vy + vz * vz);
	      double b2 = 0.5 * (bx * bx + by * by + bz * bz);
	      double pressure = et - ke - b2;
	      double ptot = et - ke;
	      pressure = pressure * gammam1;

#ifdef CHECK_NEG_PRESSURE
	      if (pressure < 0)
		{
		  std::cout
		    << __FUNCTION__
		    << ": Negative pressure"
		    << " ( " << i << " , " << j << " )" << std::endl;
		}
#endif
	      outFile << " " << pressure;
	      outFile << " " << log10 (pressure);
	      outFile << " " << log10 (b2);
	      outFile << " " << az[i][j];
	      outFile << " " << log10 (et);
	      outFile << " " << big_mesh[i][j].temperature;

	      if (big_mesh[i][j] _MASS < 0)
		{
		  std::cout
		    << __FUNCTION__ << ": "
		    << "Negative density "
		    << " (" << i << " , " << j << " ) "
		    << "dens = " << big_mesh[i][j] _MASS << std::endl;
		  MPI_Abort (MPI_COMM_WORLD, 9876);
		  exit (0);
		}
	      outFile << std::endl;
	    }
	  outFile << std::endl;
	}
      outFile.close ();

      outFile.open (blockname.c_str ());

	  for (int i = 0; i < big_x; i++)
	    {
			int j =90;
	      double gammam1 = PhysConsts::gamma - 1;
	      double rr = big_mesh[i][j] _MASS;
	      double px = big_mesh[i][j] _MOMX;
	      double py = big_mesh[i][j] _MOMY;
	      double pz = big_mesh[i][j] _MOMZ;
	      double et = big_mesh[i][j] _ENER;
	      double bx = big_mesh[i][j] _B_X;
	      double by = big_mesh[i][j] _B_Y;
	      double bz = big_mesh[i][j] _B_Z;
	      double ri = 1.0 / rr;
	      double vx = px * ri;
	      double vy = py * ri;
	      double vz = pz * ri;
	      double ke = 0.5 * rr * (vx * vx + vy * vy + vz * vz);
	      double b2 = 0.5 * (bx * bx + by * by + bz * bz);
	      double pressure = et - ke - b2;
	      double ptot = et - ke;
	      pressure = pressure * gammam1;
			outFile 
			<< " " << log10(rr)
			<< " " << log10(b2)
			<< " " << (bx)
			<< " " << (by)
			<< " " << log10(pressure)
			//<< pressure
			<< std::endl;

		 }


      outFile.close ();

//#define ZIPPER
#ifdef ZIPPER
      std::stringstream mystream;
      mystream << "/bin/gzip\t" << filename;
      std::string mystring;
      mystream >> mystring;
      std::cout << mystring.c_str () << std::endl;
      system ("/bin/gzip -f *.dat");
      //execl ("bin/ls", filename.c_str());
      //   MPI_Abort(MPI_COMM_WORLD,0);
      //  exit(0);
#endif
    }


  if (myid == 0)
    {
      std::cout << "... done" << std::endl;
    }
  delete disp;
  delete fiend;
  delete order;
  delete size;
  return 0;
}
