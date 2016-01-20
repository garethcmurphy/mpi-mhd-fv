/* $Id: ctop.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */

#include "mpi.h"
#include "out.h"
#include "PhysConsts.h"
#include "problem.h"
#include <iostream>
#include <iomanip>

int
ptoc (double *p, double *c)
{
  double rho;
  double px;
  double py, et, velx, vely, velz, velx2, vely2, velz2, ke, pressure;
  double bx, by, bz, bsquared;
  double gammag = PhysConsts::gamma;
  double gammam1 = gammag - 1;
  double gammam1i = 1 / gammam1;
  rho = p[0];
  velx = p[1];
  vely = p[2];
  velz = p[3];
  pressure = p[4];
  bx = p[5];
  by = p[6];
  bz = p[7];

  velx2 = velx * velx;
  vely2 = vely * vely;
  velz2 = velz * velz;
  ke = 0.5 * rho * (velx2 + vely2 + velz2);
  bsquared = bx * bx + by * by + bz * bz;
  c[0] = rho;
  c[1] = rho * velx;
  c[2] = rho * vely;
  c[3] = rho * velz;
  c[4] = ke + pressure * gammam1i + 0.5 * bsquared;
  c[5] = bx;
  c[6] = by;
  c[7] = bz;
  /*
     std::cout << "p " << p[0] 
     << " " << p[1] 
     << " " << p[2] 
     << " " << p[3] 
     << std::endl;
     std::cout << "c " << c[0] 
     << " " << c[1] 
     << " " << c[2] 
     << " " << c[3] 
     << std::endl;
   */

  if (c[4] <= 1e-20)
    {
      std::
	cout << "p " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] <<
	" " << p[4] << " " << p[5] << " " << p[6] << " " << p[7] << std::endl;
      std::
	cout << "c " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] <<
	" " << c[4] << " " << c[5] << " " << c[6] << " " << c[7] << std::endl;
      std::cout << "ptoc: Mein Gott!" << std::endl;
//              exit (0);
    }

  if (vely < -400)
    {
      std::cout << std::setiosflags (std::ios::fixed);
      std::
	cout << "p " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] <<
	" " << p[4] << " " << p[5] << " " << p[6] << " " << p[7] << std::endl;
      std::
	cout << "c " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] <<
	" " << c[4] << " " << c[5] << " " << c[6] << " " << c[7] << std::endl;
    }

  return 0;
}

int
ctop (double *c, double *p, double *loc)
{
  double pmin = 1e-20;
  double rho, rhoi;
  double px;
  double py, pz, et, velx, vely, velz, velx2, vely2, velz2, ke, pressure;
  double bx, by, bz, bsquared;
  double gammag = PhysConsts::gamma;
  double gammam1 = gammag - 1;
  rho = c[0];
  px = c[1];
  py = c[2];
  pz = c[3];
  et = c[4];
  bx = c[5];
  by = c[6];
  bz = c[7];
  rhoi = 1.0 / rho;
  velx = px * rhoi;
  vely = py * rhoi;
  velz = pz * rhoi;
  velx2 = velx * velx;
  vely2 = vely * vely;
  velz2 = velz * velz;
  ke = 0.5 * rho * (velx2 + vely2 + velz2);
  bsquared = 0.5 * (bx * bx + by * by + bz * bz);

  double internalenergy = et - ke - bsquared;	
  pressure = internalenergy * (gammam1);
  p[0] = rho;
  p[1] = velx;
  p[2] = vely;
  p[3] = velz;
  pressure = std::max (pressure, pmin);
  p[4] = pressure;
  p[5] = bx;
  p[6] = by;
  p[7] = bz;


#define VERBOSE_PRESSURE
#ifdef VERBOSE_PRESSURE
  if (pressure < 1e-20)
    {
      std::
	cout << "cons " << c[0] << " mv= " << c[1] << " " << c[2] << " " <<
	c[3] << " e =  " << c[4] << " b = " << c[5] << " " << c[6] << " " <<
	c[7] << std::endl;
      std::
	cout << "prim " << p[0] << " v= " << p[1] << " " << p[2] << " " <<
	p[3] << " p = " << p[4] << " b = " << p[5] << " " << p[6] << " " <<
	p[7] << std::endl;
      std::cout << __FUNCTION__ << ": pressure low" << std::endl;
      if (pressure < 0)
	{
	  std::cout << "pressure <0 " << pressure << std::endl;
	  return 1;
	}
    }
#endif
#ifdef CHECK_WEIRD_VELOCITY
  if (vely < -400)
    {
      std::cout << setiosflags (ios::fixed);
      std::
	cout << "p " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] <<
	" " << p[4] << " " << p[5] << " " << p[6] << " " << p[7] << std::endl;
      std::
	cout << "c " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] <<
	" " << c[4] << " " << c[5] << " " << c[6] << " " << c[7] << std::endl;
      std::cout << "Vel < -400" << std::endl;
    }
#endif
  return 0;
}
