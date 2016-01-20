#define MAXVAL(a,b) a>b? a:b
#include "maxspd.h"
int
maxspd (double *grid, double *maxspeed)
{


  double gammag = PhysConsts::gamma;
  double gammam1 = PhysConsts::gamma - 1;

  double ri = 0, px = 0, py = 0, pz = 0, et = 0, ke = 0;
  double rl = 0, ul = 0, vl = 0, wl, pl = 0;
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


  /* Going through every cell in the grid, 
   * work out the sound speed and take the
   * max of the absolute value of the vx and vy.
   * The sum constitutes the fastest wave speed
   * in the x or y directions and is passed back 
   * to maxspeed in the main program
   * and used to work out the time step */

  rl = grid[0];
  px = grid[1];
  py = grid[2];
  pz = 0;
  et = grid[3];

  /*
     bu = grid[5];
     bv = grid[6];
     bw = grid[7];
   */

  bu = 0;
  bv = 0;
  bw = 0;


  ri = 1.0 / rl;
  ul = px * ri;
  vl = py * ri;
  wl = pz * ri;
  vv2 = (ul * ul + vl * vl + wl * wl);
  ke = 0.5 * rl * vv2;
  bsquared = (bu * bu + bv * bv + bw * bw);
  pl = et - ke - 0.5 * bsquared;
  pl = pl * (gammam1);
  if (pl <= 0)
    {

      std::cout << __FUNCTION__ << ": pl negative, " << pl << std::endl;
      std::cout
	<< " rl = " << rl
	<< " ul = " << ul << " vl =  " << vl << " pl = " << pl << std::endl;
      std::cout << " et = " << et << std::endl;
      return (1);
    }
  bsquared = bsquared * sqrt (ri);


  calfven2 = bu * bu / sqrt (ri);
  calfven = sqrt (calfven2);
  csound2 = gammag * pl / rl;
  csound = sqrt (csound2);
  a_star2 = (csound2 + bsquared);
  term = sqrt (a_star2 * a_star2 - 4 * csound2 * calfven2);
  cfast2 = 0.5 * (a_star2 + term);
  cfast = sqrt (cfast2);

  cfast = csound;
  if (csound == 0)
    {
      std::cout << "csound = 0" << std::endl;
    }

  speed1 = fabs (ul) + cfast;
  *maxspeed = MAXVAL (speed1, *maxspeed);
//               std::cout << "speed1 " << speed1 << std::endl;

  calfven2 = bv * bv;
  term = sqrt (a_star2 * a_star2 - 4 * csound2 * calfven2);
  cfast2 = 0.5 * (a_star2 + term);
  cfast = sqrt (cfast2);
  speed2 = fabs (vl) + cfast;
  *maxspeed = MAXVAL (speed2, *maxspeed);

  calfven2 = bw * bw;
  term = sqrt (a_star2 * a_star2 - 4 * csound2 * calfven2);
  cfast2 = 0.5 * (a_star2 + term);
  cfast = sqrt (cfast2);
  speed3 = fabs (wl) + cfast;
  *maxspeed = MAXVAL (speed3, *maxspeed);


  *maxspeed =
    std::max ((std::max (fabs (ul), fabs (vl))) + csound, *maxspeed);
  if (*maxspeed < 1e-6)
    {
      std::cout << "Small speed" << std::endl;
    }
  /*
     std::cout
     << " " << rl 
     << " " << px 
     << " " << py 
     << " " << pl 
     << std::endl;
     std::cout << *maxspeed << std::endl;
     if ( *maxspeed > 700)
     exit(0);
   */


  return 0;

}
