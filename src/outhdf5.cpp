///* $Id: out.cpp,v 1.5 2006-11-16 13:48:07 gmurphy Exp $  */
/*
 *
 */

#define USE_HDF5 1
#ifdef USE_HDF5
#include "hdf5.h"
#endif


#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "tnt.h"
#include "out.h"

//using namespace std;
//using namespace TNT;

#define SCALE 1




int
outhdf5 (int numprocs, int myid, MPI_Comm new_comm, int ndims, int *dim_size,
     TNT::Array2D < unk > small_mesh, int *arr_size, int step)
{
  const int ne = NE;
  const int small_x = arr_size[0];
  const int small_y = arr_size[1];
  TNT::Array2D < unk > spare_mesh (200, 200);
  TNT::Array2D < double >small_arr (small_x, small_y);
  //int big_x = 9;
  //int big_y = 9;
  const int big_x = (small_x - 2 * NGC) * (dim_size[0]);
  const int big_y = (small_y - 2 * NGC) * (dim_size[1]);
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
  MPI_Aint displacements[3] = { 0, (int) test.array - (int) &test.temperature, sizeof (test) };	/* guessing... */
  int blks[3] = { 1, ne, 1 };
  MPI_Datatype types_str[3] = { MPI_INT, MPI_DOUBLE, MPI_UB };


  MPI_Type_struct (3, blks, displacements, types_str, &mycolumn);
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

  int sizeofdouble = 0;
  int sizeofunk = 0;
  MPI_Type_extent (MPI_DOUBLE, (MPI_Aint *) & sizeofdouble);
  MPI_Type_extent (mycolumn, (MPI_Aint *) & sizeofunk);

  // Define a datatype with a stride 
//  MPI_Type_vector (small_x, 1, 2, MPI_DOUBLE, &small_row);
//  MPI_Type_hvector (small_y, 1, big_x*sizeofdouble, MPI_DOUBLE, &small_block);

  MPI_Type_hvector ((small_x - 2 * NGC), (small_y - 2 * NGC),
		    big_y * sizeofunk, mycolumn, &small_row);
  MPI_Type_commit (&small_row);

  MPI_Aint disps[2] = { 0, (small_x - 2 * NGC) * sizeof (unk) };	/* guessing... */
  int blocklengths[2] = { 1, 1 };
  MPI_Datatype types[2] = { small_row, MPI_UB };

  MPI_Type_struct (2, blocklengths, disps, types, &big_type);
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

#define GNUPLOT
#ifdef GNUPLOT

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
#endif

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



#ifdef USE_HDF5
#define RANK 2
  hid_t file, dataset;     /* file and dataset handles */
  hid_t datatype, dataspace;  /* handles */
  hsize_t dimsf[2];     /* dataset dimensions */
  herr_t status;
  TNT::Array2D <double> data( big_x,big_y);
  //double data[ (const int)big_x][(const int)big_y];
  //double data[200][200];
  double *bufferArray= (double * ) malloc( sizeof(double ) *big_x*big_y);
#ifdef ZINGALE
  double ** data= (double **) malloc( sizeof(double *) *big_y);

  for (int i = 0; i < big_y; i++) {
	      data[i] = bufferArray + i*big_x;
			  }
#endif

  char *names[] =
    { "Density", "Velx", "Vely", "Velz", "Energy", "Bx", "By", "Bz" };
  int ll = 0;
  std::stringstream hdf5_stream_filename;
  std::string hdf5_filename;

      std::string str_file_tag2;
      stream_tag.clear ();
      stream_tag.width (4);
      stream_tag.fill ('0');
      stream_tag << step;
      stream_tag >> str_file_tag2;
  hdf5_stream_filename << outputdir << "hdf5_"  << str_file_tag2 << ".h5";
  hdf5_stream_filename >> hdf5_filename;
	std::cout << hdf5_filename.c_str() << std::endl;
  file = H5Fcreate (hdf5_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
          H5P_DEFAULT);
#define WRITEDATA
#ifdef WRITEDATA
  for (ll = 0; ll < 8; ll++)
    { 
      dimsf[0] = big_x;
      dimsf[1] = big_y;
      dataspace = H5Screate_simple (RANK, dimsf, NULL);
      /*
       * Define datatype for the data in the file.
       * We will store little endian DOUBLE numbers.
       */
      datatype = H5Tcopy (H5T_NATIVE_DOUBLE);
      status = H5Tset_order (datatype, H5T_ORDER_LE);
      /*
       * Create a new dataset within the file using defined dataspace and
       * datatype and default dataset creation properties.
       */
      dataset = H5Dcreate (file, names[ll], datatype, dataspace, H5P_DEFAULT);

      for (int jj = 0; jj < big_y; jj++)
   {
     for (int ii = 0; ii < big_x; ii++)
	  {
      // data[ii][jj] = ii*jj*1.0;
       data[ii][jj] = big_mesh[ii][jj].array[ll];
		 //std::cout << ii << " " << jj << std::endl;
	  }
   }
      /*
       * Write the data to the dataset using default transfer properties.
       */
      //status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
		// std::cout << "just before write "  << std::endl;
      status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(data[0][0]));

      /*
       * Close/release resources.
       */
      H5Sclose (dataspace);
      H5Tclose (datatype);
      H5Dclose (dataset);
    }


#endif

  H5Fclose (file);

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
