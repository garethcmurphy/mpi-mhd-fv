/* $Id: output.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "output.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


//#define USE_HDF5 1
#ifdef USE_HDF5
#include "hdf5.h"
#endif


#define H5FILE_NAME        "SDS.h5"
#define RANK   2

#define DATASETNAME "IntArray"


int
output (Array3D < zone > grid, Array3D < zone > fx, Array3D < zone > fy,
	int time, char *filename)
{

#ifdef USE_HDF5
  hid_t file, dataset;		/* file and dataset handles */
  hid_t datatype, dataspace;	/* handles */
  hsize_t dimsf[2];		/* dataset dimensions */
  herr_t status;
  double data[nx][ny];
  char *names[] =
    { "Density", "Velx", "Vely", "Velz", "Energy", "Bx", "By", "Bz" };
  int ll = 0;
  stringstream hdf5_stream_filename;
  string hdf5_filename;
#endif


  ofstream fout;
  ofstream gout;
  double gammam1 = gammag - 1;
  double rl, ri;
  double px;
  double py;
  double pz;
  double pl;
  double bx;
  double by;
  double bz;
  double bsquared;
  double et, ul, vl, wl, ke, al;
  int ii = 0;
  int jj = 0;
  int kk = 0;
  char outputdir[50] = "output/";
//      char            filename[50] = "out_2d_";
  stringstream s;
  stringstream stream_filename;
  stringstream stream_temp_b;
  string str_file_tag;
  string str_output_filename;
  string str_input_filename;


  double ki = 24296.3696;
  double mp = 1.67262158;
  double mpi = 1.0 / mp;
  double nt = 0;
  double nt2 = 0;
  double temperature = 0;

  s.clear ();
  s.width (5);
  s.fill ('0');
  s << time;
  s >> str_file_tag;
  stream_filename.clear ();
  stream_filename << outputdir << filename << str_file_tag;
  stream_filename >> str_input_filename;

#ifdef USE_HDF5
  hdf5_stream_filename << outputdir << "hdf5_" << filename << str_file_tag <<
    ".h5";
  hdf5_stream_filename >> hdf5_filename;
  file =
    H5Fcreate (hdf5_filename.c_str (), H5F_ACC_TRUNC, H5P_DEFAULT,
	       H5P_DEFAULT);

  for (ll = 0; ll < ne; ll++)
    {
      dimsf[0] = nx;
      dimsf[1] = ny;
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

      for (jj = 0; jj < ny; jj++)
	{
	  for (ii = 0; ii < nx; ii++)
	    data[ii][jj] = grid[ii][jj][kk].array[ll];
	}
      /*
       * Write the data to the dataset using default transfer properties.
       */
      status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
			 H5P_DEFAULT, data);

      /*
       * Close/release resources.
       */
      H5Sclose (dataspace);
      H5Tclose (datatype);
      H5Dclose (dataset);
    }


  dimsf[0] = nx;
  dimsf[1] = ny;
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
  dataset = H5Dcreate (file, "Pressure", datatype, dataspace, H5P_DEFAULT);

  for (jj = 0; jj < ny; jj++)
    {
      for (ii = 0; ii < nx; ii++)
	{

	  rl = grid[ii][jj][kk] _MASS;
	  px = grid[ii][jj][kk] _MOMX;
	  py = grid[ii][jj][kk] _MOMY;
	  pz = grid[ii][jj][kk] _MOMZ;
	  et = grid[ii][jj][kk] _ENER;
	  bx = grid[ii][jj][kk] _B_X;
	  by = grid[ii][jj][kk] _B_Y;
	  bz = grid[ii][jj][kk] _B_Z;
	  ri = 1.0 / rl;
	  ul = px * ri;
	  vl = py * ri;
	  wl = pz * ri;
	  ke = 0.5 * rl * (ul * ul + vl * vl + wl * wl);
	  bsquared = bx * bx + by * by + bz * bz;
	  pl = et - ke - 0.5 * bsquared;
	  pl = pl * gammam1;
	  al = sqrt (gammag * pl * ri);
	  nt = 2 * mpi * rl;
	  nt2 = nt * nt;
	  temperature = ki * pl / nt;
	  temperature = log10 (temperature);


	  data[ii][jj] = pl;
	}
    }
  /*
   * Write the data to the dataset using default transfer properties.
   */
  status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, data);

  /*
   * Close/release resources.
   */
  H5Sclose (dataspace);
  H5Tclose (datatype);
  H5Dclose (dataset);


  H5Fclose (file);
//#else 


  fout.open (str_input_filename.c_str ());
  if (!fout)
    {
      cerr << "unable to open file " << endl;
    }

  jj = 0;
  for (ii = 0; ii < nx; ii++)
    {
//       for (jj = 0; jj < ny; jj++)
      {
	rl = grid[ii][jj][kk] _MASS;
	px = grid[ii][jj][kk] _MOMX;
	py = grid[ii][jj][kk] _MOMY;
	pz = grid[ii][jj][kk] _MOMZ;
	et = grid[ii][jj][kk] _ENER;
	bx = grid[ii][jj][kk] _B_X;
	by = grid[ii][jj][kk] _B_Y;
	bz = grid[ii][jj][kk] _B_Z;
	ri = 1.0 / rl;
	ul = px * ri;
	vl = py * ri;
	wl = pz * ri;
	ke = 0.5 * rl * (ul * ul + vl * vl + wl * wl);
	bsquared = bx * bx + by * by + bz * bz;
	pl = et - ke - 0.5 * bsquared;
	pl = pl * gammam1;
	al = sqrt (gammag * pl * ri);

#ifdef DEBUG_BC
	if (ii == 2 && jj == 2 && px != 0)
	  {
	    cout << px << endl;
	    cout << ul << endl;
	    cout << "wtf?" << endl;
	  }
#endif /* DEBUG_BC */

	fout
	  << " " << rl
	  << " " << ul
	  << " " << vl
	  << " " << wl
	  << " " << et
	  << " " << bx
	  << " " << by << " " << bz << " " << pl << " " << al << endl;
      }
    }
  fout.close ();

  gout.open ("gm.general");
  gout << "file = /home/gmurphy/mhdvanleer-0.0.1/" << str_input_filename <<
    endl;
  gout << "grid = " << nx << " x " << ny << endl;
  gout << "format = ascii" << endl;
  gout << "interleaving = field" << endl;
  gout << "majority = row" << endl;
  gout << "field = V_sound, E_tot, Rho, Vel_X, Vel_Y, Pressure" << endl;
  gout << "structure = scalar, scalar, scalar, scalar, scalar, scalar" <<
    endl;
  gout << "type = double, double, double, double, double, double" << endl;
  gout <<
    "dependency = positions, positions, positions, positions, positions, positions"
    << endl;
  gout << "positions = regular, regular, 0, 1, 0, 1" << endl;
  gout << "" << endl;
  gout << "end" << endl;
  gout.close ();

#endif /* HDF5 or not  */
  return 0;
}




int
write_data_to_hdf5_file (int nx, int ny, double **data, hid_t file)
{

  hid_t dataset;		/* file and dataset handles */
  hid_t datatype, dataspace;	/* handles */
  hsize_t dimsf[2];		/* dataset dimensions */
  herr_t status;


  dimsf[0] = nx;
  dimsf[1] = ny;
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
  dataset = H5Dcreate (file, "Temperature", datatype, dataspace, H5P_DEFAULT);

  /*
   * Write the data to the dataset using default transfer properties.
   */
  status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, data);

  /*
   * Close/release resources.
   */
  H5Sclose (dataspace);
  H5Tclose (datatype);
  H5Dclose (dataset);

  return 0;
}
