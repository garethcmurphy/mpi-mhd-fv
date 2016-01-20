/* $Id: log.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "log.h"
int
log (Array3D < zone > grid, Array3D < zone > fx, Array3D < zone > fy,
     double timestep, double maximumspeed, double time, double del)
{

  int ii = 0;
  int jj = 0;
  int kk = 0;

  ofstream logfile;

  double px, py, et, ri, rl, ul, vl, ke, pl, al;
  double gammam1 = gammag - 1;



  logfile.open ("logfile.txt", ios::app);
  logfile << "Timestep=" << setiosflags (ios::fixed) << timestep;
  logfile << " Maxspeed= " << setiosflags (ios::fixed) << maximumspeed;
  logfile << " Time=" << setiosflags (ios::fixed) << time;
  logfile << " cfl= " << setiosflags (ios::scientific) << del;
  logfile << endl;
  //     for (jj = 1; jj < ny; jj++)
//     {
  for (ii = 1; ii < nx - 1; ii++)
    {
      logfile
	<< ii
	<< " "
	<< jj
	<< " rho\t"
	<< setiosflags (ios::fixed)
	<< setw (9)
	<< setprecision (7)
	<< grid[ii][jj][kk] _MASS
	<< "\t" << del
	<< "\t" << fx[ii + 1][jj][kk] _MASS
	<< "\t" << fx[ii][jj][kk] _MASS << endl;
      logfile
	<< ii
	<< " " << jj
	<< " vx\t"
	<< setiosflags (ios::fixed)
	<< setw (5)
	<< setprecision (7)
	<< grid[ii][jj][kk] _MOMX
	<< "\t" << del
	<< "\t" << fx[ii + 1][jj][kk] _MOMX
	<< "\t" << fx[ii][jj][kk] _MOMX << endl;
      logfile
	<< ii
	<< " "
	<< jj
	<< " vy\t"
	<< setiosflags (ios::fixed)
	<< setw (5)
	<< setprecision (7)
	<< grid[ii][jj][kk] _MOMY
	<< "\t" << del
	<< "\t" << fx[ii + 1][jj][kk] _MOMY
	<< "\t" << fx[ii][jj][kk] _MOMY << endl;

      logfile
	<< ii
	<< " "
	<< jj
	<< " ene\t" << setiosflags (ios::
				    fixed) << setw (9) <<
	setprecision (7) << grid[ii][jj][kk] _ENER
	<< "\t" << del
	<< "\t" << fx[ii + 1][jj][kk] _ENER
	<< "\t" << fx[ii][jj][kk] _ENER << endl;

      rl = grid[ii][jj][kk] _MASS;
      px = grid[ii][jj][kk] _MOMX;
      py = grid[ii][jj][kk] _MOMY;
      et = grid[ii][jj][kk] _ENER;
      ri = 1.0 / rl;
      ul = px * ri;
      vl = py * ri;
      ke = 0.5 * rl * (ul * ul + vl * vl);
      pl = et - ke;
      pl = pl * gammam1;
      al = sqrt (gammag * pl * ri);

      logfile
	<< ii
	<< " "
	<< jj
	<< " p\t" << setiosflags (ios::
				  fixed) << setw (9) <<
	setprecision (7) << pl << endl;
      logfile
	<< ii
	<< " "
	<< jj
	<< " a\t" << setiosflags (ios::
				  fixed) << setw (9) <<
	setprecision (7) << al << endl;

      logfile << endl;

    }
//     }

  logfile.close ();

  return 0;
}
