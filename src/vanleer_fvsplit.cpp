/* $Id: vanleer_fvsplit.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "vanleer_fvsplit.h"
int
vanleer_flux_vector_split (double *leftstate, double *rightstate, double *fp,
			   double *fn, int iii, int jjj, int idir)
{
  int hh = 0;

  double gammam1 = gammag - 1;
  double gammai = 1 / gammag;
  double rhoi, px, py, et, ke, p;
  double rho, velx, vely, pres, vsnd, mach;
  double velx2;
  double vely2;
  double massfluxp;
  double massfluxn;

  /* Determine the Mach number */

  rho = leftstate[0];
  px = leftstate[1];
  py = leftstate[2];
  et = leftstate[3];


  rhoi = 1.0 / rho;
  velx = px * rhoi;
  vely = py * rhoi;
  velx2 = velx * velx;
  vely2 = vely * vely;
  ke = 0.5 * rho * (velx2 + vely2);
  pres = et - ke;
  pres = pres * (gammam1);
  vsnd = gammag;
  vsnd = vsnd * rhoi * pres;
  vsnd = sqrt (vsnd);

  if (idir == 1)
    {
      mach = velx / vsnd;
      massfluxp = rho * pow ((velx + vsnd), 2) / (4 * vsnd);
      massfluxn = (-rho) * pow ((velx - vsnd), 2) / (4 * vsnd);
      /* If Mach <= -1 */
      if (mach <= -1)
	{
	  /* Find f+ and f- */
	  fp[0] = 0;
	  fp[1] = 0;
	  fp[2] = 0;
	  fp[3] = 0;

	  fn[0] = rho * velx;
	  fn[1] = rho * velx2 + pres;
	  fn[2] = rho * velx * vely;
	  fn[3] = velx * (et + p);
	}
      /* If Mach > 1 */
      else if (mach >= 1)
	{
	  /* Find f+ and f- */
	  fp[0] = rho * velx;
	  fp[1] = rho * velx2 + pres;
	  fp[2] = rho * velx * vely;
	  fp[3] = velx * (et + p);

	  fn[0] = 0;
	  fn[1] = 0;
	  fn[2] = 0;
	  fn[3] = 0;
	}

      else
	{
	  /* If -1< Mach < 1 */
	  /* Find f+ and f- */
	  fp[0] = massfluxp;
	  fp[1] = massfluxp * (gammam1 * velx + 2.0 * vsnd) * gammai;
	  fp[2] = massfluxp * vely;
	  fp[3] =
	    massfluxp * (vely2 * 0.5 +
			 pow ((gammam1 * velx + 2.0 * vsnd),
			      2) * 0.5 / (gammag * gammag - 1));

	  fn[0] = massfluxn;
	  fn[1] = massfluxn * (gammam1 * velx - 2.0 * vsnd) * gammai;
	  fn[2] = massfluxn * vely;
	  fn[3] =
	    massfluxn * (vely2 * 0.5 -
			 pow ((gammam1 * velx - 2.0 * vsnd),
			      2) * 0.5 / (gammag * gammag - 1));
	}
      /* Sum them and return them in one flux */
      /* Not too sure about this ... */
    }


  else if (idir == 2)
    {
      mach = vely / vsnd;
      massfluxp = rho * pow ((vely + vsnd), 2) / (4 * vsnd);
      massfluxn = (-rho) * pow ((vely - vsnd), 2) / (4 * vsnd);
      /* If Mach <= -1 */
      if (mach <= -1)
	{
	  /* Find f+ and f- */
	  fp[0] = 0;
	  fp[1] = 0;
	  fp[2] = 0;
	  fp[3] = 0;

	  fn[0] = rho * vely;
	  fn[1] = rho * velx * vely;
	  fn[2] = rho * vely2 + pres;
	  fn[3] = vely * (et + p);
	}
      /* If Mach > 1 */
      else if (mach >= 1)
	{
	  /* Find f+ and f- */
	  fp[0] = rho * vely;
	  fp[1] = rho * vely * velx;
	  fp[2] = rho * vely2 + pres;
	  fp[3] = vely * (et + p);

	  fn[0] = 0;
	  fn[1] = 0;
	  fn[2] = 0;
	  fn[3] = 0;
	}

      else
	{
	  /* If -1< Mach < 1 */
	  /* Find f+ and f- */
	  fp[0] = massfluxp;
	  fp[1] = massfluxp * velx;
	  fp[2] = massfluxp * (gammam1 * vely + 2 * vsnd) * gammai;
	  fp[3] =
	    massfluxp * (velx2 * 0.5 +
			 pow ((gammam1 * vely + 2 * vsnd),
			      2) * 0.5 / (gammag * gammag - 1));

	  fn[0] = massfluxn;
	  fn[1] = massfluxn * velx;
	  fn[2] = massfluxn * (gammam1 * vely - 2 * vsnd) * gammai;
	  fn[3] =
	    massfluxn * (velx2 * 0.5 -
			 pow ((gammam1 * vely - 2 * vsnd),
			      2) * 0.5 / (gammag * gammag - 1));
	}
      /* Sum them and return them in one flux */
      /* Not too sure about this ... */
    }

//      cout << "vanl: " << fp[0] << endl;
  return 0;

}
