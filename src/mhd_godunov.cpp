/* $Id: mhd_godunov.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
/* $Id: mhd_godunov.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $: mhd_godunov.cpp,v 1.2 2006-10-30 15:12:38 gmurphy Exp $ */
#include "mhd_godunov.h"
#undef DEBUG
int
godunov_mhd (Array3D < zone > gr,
	     Array3D < zone > fl,
	     Array3D < zone > flux2,
	     int height, int width, int jj, double *time)
{
  ofstream outFile;
  char outfilename[50] = "output/godunov.log";
  int ii = 0;
  int kk = 0;
  int hh = 0;
  int rc = 0;
  int idir = 0;
  double leftstate[8];
  double rightstate[8];
  double k = 0;
  double max_speed = 0;
  double temp = 0;
  double delta_t = 0;
  double cfl = 0.8;
//  k=0.25/4.44;
//  k=0.3/4.44;

  /*
   * Determine fastest wave speed, time step and k
   * */
  max_speed = 0;
  for (ii = 1; ii <= height - 1; ii++)
    {
      leftstate[hh] = gr[ii - 1][jj - 1][kk].array[hh];
      rightstate[hh] = gr[ii][jj - 1][kk].array[hh];
      temp = max_speed;
      rc =
	falle_mhd (leftstate, rightstate, fl[ii - 1][jj - 1][kk].array, jj,
		   &max_speed, idir);
      if (temp > max_speed)
	{
	  max_speed = temp;
	}
    }

  k = cfl / max_speed;
  delta_t = cfl / (max_speed * 512);
  *time = *time + delta_t;
  cout
    << " Timestep= " << jj
    << " Time= " << *time
    << " Max speed= " << max_speed << " CFL= " << k << endl;

  for (hh = 1; hh < height - 1; hh++)
    {
      gr[ii][jj][kk].array[hh] =
	gr[ii][jj - 1][kk].array[hh] - k * (fl[ii][jj - 1][kk].array[hh] -
					    fl[ii - 1][jj - 1][kk].array[hh]);
    }

#ifdef DEBUG
  if (1)
    {
      outFile.open (outfilename, ofstream::out | ofstream::app);

      if (!outFile)
	{
	  cerr << "mhd:Unable to open file" << endl;
	}

      outFile << "jj " << jj << endl;

      for (ii = 1; ii < height - 1; ii++)
	{
	  if (fl[ii][jj - 1].mass != fl[ii - 1][jj - 1].mass)
	    {
	      outFile << ii << " rho " << gr[ii][jj].mass << " "
		<< gr[ii][jj - 1].mass << " "
		<< fl[ii][jj - 1].mass << " "
		<< fl[ii - 1][jj - 1].mass << endl;
	      outFile << ii << " px  " << gr[ii][jj].px << " "
		<< gr[ii][jj - 1].px << " "
		<< fl[ii][jj - 1].px << " " << fl[ii - 1][jj - 1].px << endl;
	      outFile << ii << " py  " << gr[ii][jj].py << " "
		<< gr[ii][jj - 1].py << " "
		<< fl[ii][jj - 1].py << " " << fl[ii - 1][jj - 1].py << endl;
	      outFile << ii << " pz  " << gr[ii][jj].pz << " "
		<< gr[ii][jj - 1].pz << " "
		<< fl[ii][jj - 1].pz << " " << fl[ii - 1][jj - 1].pz << endl;
	      outFile << ii << " bx  " << gr[ii][jj].bx << " "
		<< gr[ii][jj - 1].bx << " "
		<< fl[ii][jj - 1].bx << " " << fl[ii - 1][jj - 1].bx << endl;
	      outFile << ii << " by  " << gr[ii][jj].by << " "
		<< gr[ii][jj - 1].by << " "
		<< fl[ii][jj - 1].by << " " << fl[ii - 1][jj - 1].by << endl;
	      outFile << ii << " bz  " << gr[ii][jj].bz << " "
		<< gr[ii][jj - 1].bz << " "
		<< fl[ii][jj - 1].bz << " " << fl[ii - 1][jj - 1].bz << endl;
	      outFile << ii << " ene " << gr[ii][jj].energy << " "
		<< gr[ii][jj - 1].energy << " "
		<< fl[ii][jj - 1].energy << " "
		<< fl[ii - 1][jj - 1].energy << endl;
	      outFile << endl;
	    }
	}
      if (!outFile)
	{
	  cerr << "mhd:Unable to open file" << endl;
	}
    }
#endif /* DEBUG */
  return 0;
}
