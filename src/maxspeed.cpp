/* $Id: maxspeed.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "mpi.h"
#include "maxspeed.h"
int
maxspeed (TNT::Array2D < unk > mesh, double *grid, double *maxspeed)
{
  int ii = 0;
  int jj = 0;
  int kk = 0;


  double gammag = 5.0 / 3.0;
  double gammam1 = gammag - 1;

  double ri = 0, px = 0, py = 0, pz = 0, et = 0, ke = 0;
  double rl = 0, ul = 0, vl = 0, wl, pl = 0, al = 0, speed = 0;
  double bu, bv, bw, bsquared;
  double calfven = 0;
  double calfven2 = 0;
  double cfast = 0;
  double cfast2 = 0;
  double csound2 = 0;
  double csound = 0;
  double term = 0;
  double a_star2 = 0;
  double vv2 = 0;
  double speed1 = 0;
  double speed2 = 0;
  double speed3 = 0;


  /* Going through every cell in the grid, work out the sound speed and take the
   * max of the absolute value of the vx and vy. The sum constitutes the fastest wave speed
   * in the x or y directions and is passed back to maxspeed in the main program
   * and used to work out the time step */

  int nx = mesh.dim1 ();
  int ny = mesh.dim2 ();

  int maxi = 0;
  int maxj = 0;

  for (ii = 2; ii < nx - 2; ii++)
    {
      for (jj = 2; jj < ny - 2; jj++)
	{

	  rl = mesh[ii][jj].array[0];
	  px = mesh[ii][jj].array[1];
	  py = mesh[ii][jj].array[2];
	  pz = mesh[ii][jj].array[3];
	  et = mesh[ii][jj].array[4];

	  bu = mesh[ii][jj].array[5];
	  bv = mesh[ii][jj].array[6];
	  bw = mesh[ii][jj].array[7];



	  ri = 1.0 / rl;
	  ul = px * ri;
	  vl = py * ri;
	  wl = pz * ri;
	  vv2 = (ul * ul + vl * vl + wl * wl);
	  ke = 0.5 * rl * vv2;
	  bsquared = (bu * bu + bv * bv + bw * bw);
	  pl = et - ke - 0.5 * bsquared;
	  pl = pl * (gammam1);
	  bsquared = bsquared * ri;


	  calfven2 = bu * bu * ri;
	  calfven = sqrt (calfven2);
	  csound2 = gammag * pl / rl;
	  csound = sqrt (csound2);
	  a_star2 = (csound2 + bsquared);
	  term = sqrt (a_star2 * a_star2 - 4 * csound2 * calfven2);
	  cfast2 = 0.5 * (a_star2 + term);
	  cfast = sqrt (cfast2);


	  speed = std::max (fabs (ul), fabs (vl)) + cfast;
	  speed1 = fabs (ul) + cfast;
//               speed1 = fabs(ul)+csound+calfven;

	  if (speed1 > *maxspeed)
	    {
	      maxi = ii;
	      maxj = jj;
	    }

	  *maxspeed = std::max (speed1, *maxspeed);

	  calfven2 = bv * bv * ri;
	  term = sqrt (a_star2 * a_star2 - 4 * csound2 * calfven2);
	  cfast2 = 0.5 * (a_star2 + term);
	  cfast = sqrt (cfast2);
	  speed2 = fabs (vl) + cfast;
	  if (speed2 > *maxspeed)
	    {
	      maxi = ii;
	      maxj = jj;
	    }
	  *maxspeed = std::max (speed2, *maxspeed);

	  calfven2 = bw * bw * ri;
	  term = sqrt (a_star2 * a_star2 - 4 * csound2 * calfven2);
	  cfast2 = 0.5 * (a_star2 + term);
	  cfast = sqrt (cfast2);
	  speed3 = fabs (wl) + cfast;
	  *maxspeed = std::max (speed3, *maxspeed);

	  if (isinf (*maxspeed))
	    {
	      std::cout
			<< std::endl
		<< __FUNCTION__ << " : "
		<< " ( " << ii << " , " << jj << " ) " << std::endl;
	      MPI_Abort (MPI_COMM_WORLD, 3456789);
	      exit (0);

	    }

	}
    }

//          *maxspeed = max (speed, *maxspeed);

  /*
     cout << ii
     << " " << jj
     << " " << rl 
     << " " << px 
     << " " << py 
     << " " << pl 
     << endl;
     cout << *maxspeed << endl;
     if ( *maxspeed > 700)
     exit(0);
   */

#define WHEREISMAXSPEED
#ifdef WHEREISMAXSPEED
  int myid;
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  std::cout
    << "\nProc " << myid << " Max speed location (" << maxi
    << "," << maxj << ")" << " =  " << *maxspeed << std::endl;
#endif

  return 0;

}
