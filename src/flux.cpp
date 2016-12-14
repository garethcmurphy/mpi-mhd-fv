/* $Id: flux.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "tnt.h"
#define ROE
#undef VANLEER
#define SECOND_ORDER_SPACE
//#undef SECOND_ORDER_SPACE
int ctop (double *c, double *p);
int ptoc (double *p, double *c);
int
flux (TNT::Array2D < unk > oldgrid,
      TNT::Array2D < flux > fx,
      TNT::Array2D < flux > fy,
      double *InterfaceFlux,
      double *ResolvedState,
      double dtodx, int ii, int jj, int timestep, int idir, int second)
{
  int d[2];
  int kk = 0;
  int hh = 0;
  int rc;

  double slope1 = 0.0;
  double slope2 = 0.0;
  double left = 0.0;
  double mid = 0.0;
  double right = 0.0;

  double *leftstate;
  double *rightstate;
  double *flux;
  double *fp;
  double *fn;


  double pll[ne];
  double plm[ne];
  double plr[ne];
  double prr[ne];
  double leftprim[ne];
  double rightprim[ne];
  double xx[ne];
  double yy[ne];
  double unused;


  leftstate = new double[ne];
  rightstate = new double[ne];
  flux = new double[ne];
  fp = new double[ne];
  fn = new double[ne];


  if (idir == 1)
    {
      d[0] = 1;
      d[1] = 0;
    }
  else if (idir == 2)
    {
      d[0] = 0;
      d[1] = 1;
    }

  /* Determine the x fluxes */
  for (hh = 0; hh < ne; hh++)
    {
      leftstate[hh] = oldgrid[ii - d[0]][jj - d[1]][kk].array[hh];
      rightstate[hh] = oldgrid[ii][jj][kk].array[hh];
//               leftstate[hh]  = fabs(leftstate[hh]) > 1e-10?  leftstate[hh] :0; 
//               rightstate[hh]  = fabs(rightstate[hh]) > 1e-10?  rightstate[hh] :0; 
//               cout << "flux:" << leftstate[hh] << " " << rightstate[hh] << endl;
    }

  rc = ctop (leftstate, leftprim);
  rc = ctop (rightstate, rightprim);
//     Need a 2nd order in space correction and a non-linear flux limiter here

#ifdef SECOND_ORDER_SPACE
  if (second == 1)
//  if (1)
    {

      rc = ctop (oldgrid[ii - 2 * d[0]][jj - 2 * d[1]][kk].array, pll);
      /*
         if (rc == 1)
         {
         cout << "Location: " << ii << ", " << jj << endl;
         exit (0);
         }
       */
      pll[4] = max (pmin, pll[4]);
      rc = ctop (oldgrid[ii - d[0]][jj - d[1]][kk].array, plm);
      /*
         if (rc == 1)
         {
         cout << "Location: " << ii << ", " << jj << endl;
         exit (0);
         }
       */
      plm[4] = max (pmin, plm[4]);
      rc = ctop (oldgrid[ii][jj][kk].array, plr);
      /*
         if (rc == 1)
         {
         cout << "Location: " << ii << ", " << jj << endl;
         exit (0);
         }
       */
      plr[4] = max (pmin, plr[4]);
      /*
         rc = ctop (oldgrid[ii + d[0]][jj + d[1]][kk].array, prr);
         if (rc == 1)
         {
         cout << "Location: " << ii << ", " << jj << endl;
         exit (0);
         }
       */
      prr[4] = max (pmin, prr[4]);

      for (hh = 0; hh < ne; hh++)
	{
//               left  = oldgrid[ii - 2 * d[0]][jj - 2 * d[1]][kk].array[hh];
//               mid   = oldgrid[ii -     d[0]][jj - d[1]][kk].array[hh];
//               right = oldgrid[ii           ][jj][kk].array[hh];
	  left = pll[hh];
	  mid = plm[hh];
	  right = plr[hh];
	  slope1 = vanleer ((mid - left), (right - mid));

//               left  = oldgrid[ii - d[0]][jj - d[1]][kk].array[hh];
//               mid   = oldgrid[ii       ][jj       ][kk].array[hh];
//               right = oldgrid[ii + d[0]][jj + d[1]][kk].array[hh];
	  left = plm[hh];
	  mid = plr[hh];
	  right = prr[hh];
	  slope2 = vanleer ((mid - left), (right - mid));

	  leftprim[hh] = plm[hh] + 0.5 * delta_x * slope1;
	  rightprim[hh] = plr[hh] - 0.5 * delta_x * slope2;
	}

      // Flatten slopes here?

      leftprim[4] = max (pmin, leftprim[4]);
      rc = ptoc (leftprim, leftstate);
/*		if (rc ==1) {
			cout << __FUNCTION__ << ": Exiting on neg pressure." << endl;
			exit(0);
		}
		*/
      rightprim[4] = max (pmin, rightprim[4]);
      rc = ptoc (rightprim, rightstate);
      /*
         if (rc ==1) {
         cout << __FUNCTION__ << ": Exiting on neg pressure." << endl;
         exit(0);
         }
       */
      if (leftprim[4] < 1e-10)
	{
	  cout << "leftprim[4] < 1e-10 " << leftprim[4]
	    << " " << oldgrid[ii - d[0]][jj - d[1]][kk].array[hh]
	    << " " << ii << " " << jj << endl;
	  if (leftprim[4] < 0)
	    {
	      cout << " Exiting on pressure, " << leftprim[4] << " < 0 " <<
		endl;
	      exit (0);
	    }
	}
      if (rightprim[4] < 1e-10)
	{
	  cout << "rightprim[4] < 1e-10 " << rightprim[4]
	    << " " << oldgrid[ii][jj][kk].array[hh]
	    << " " << ii << " " << jj << endl;
	  if (rightprim[4] < 0)
	    {
	      cout << " Exiting on pressure, " << rightprim[4] << " < 0 " <<
		endl;
	      exit (0);
	    }
	}


    }
#endif


#ifdef ROE
  if (leftstate[0] < 0 || rightstate[0] < 0)
    {
      cout << "flux:negative density in roe solver" << endl;
      cout << "Location = " << ii << "," << jj << endl;
      exit (0);
    }

  // Rotate fluxes here
  double lstate[8];
  double rstate[8];
  double iflux[8];
  double Res_state[8];

  if (idir == 1)
    {
      for (hh = 0; hh < ne; hh++)
	{
	  lstate[hh] = leftstate[hh];
	  rstate[hh] = rightstate[hh];
	}
    }
  else if (idir == 2)
    {
      rstate[0] = rightstate[0];
      rstate[1] = rightstate[2];
      rstate[2] = rightstate[3];
      rstate[3] = rightstate[1];
      rstate[4] = rightstate[4];
      rstate[5] = rightstate[6];
      rstate[6] = rightstate[7];
      rstate[7] = rightstate[5];

      lstate[0] = leftstate[0];
      lstate[1] = leftstate[2];
      lstate[2] = leftstate[3];
      lstate[3] = leftstate[1];
      lstate[4] = leftstate[4];
      lstate[5] = leftstate[6];
      lstate[6] = leftstate[7];
      lstate[7] = leftstate[5];
    }


//rc = riemann (leftstate, rightstate, InterfaceFlux, ResolvedState, timestep, &unused, idir);
//  rc = riemann (lstate, rstate, iflux, Res_state, timestep, &unused, 1);
//  rc = hlld (leftprim, rightprim, iflux, Res_state, timestep, &unused, 1);
  rc = riemann (lstate, rstate, iflux, Res_state, timestep, &unused, 1);

  if (rc == 1)
    {
      cout << __FUNCTION__ << ": location " << ii << "," << jj << endl;
      exit (0);
    }

  if (idir == 1)
    {
      for (hh = 0; hh < ne; hh++)
	{
	  InterfaceFlux[hh] = iflux[hh];
	  ResolvedState[hh] = Res_state[hh];
	}
    }
  else if (idir == 2)
    {
      InterfaceFlux[0] = iflux[0];
      InterfaceFlux[2] = iflux[1];
      InterfaceFlux[3] = iflux[2];
      InterfaceFlux[1] = iflux[3];
      InterfaceFlux[4] = iflux[4];
      InterfaceFlux[6] = iflux[5];
      InterfaceFlux[7] = iflux[6];
      InterfaceFlux[5] = iflux[7];

      ResolvedState[0] = Res_state[0];
      ResolvedState[2] = Res_state[1];
      ResolvedState[3] = Res_state[2];
      ResolvedState[1] = Res_state[3];
      ResolvedState[4] = Res_state[4];
      ResolvedState[6] = Res_state[5];
      ResolvedState[7] = Res_state[6];
      ResolvedState[5] = Res_state[7];
    }

  // Unrotate fluxes


  // Artificial Viscosity goes in here  ----
  //
  // LAPIDUS method divv()*q * Conserved Differences* 
  // Second-order divv either at x-face or y-face
  //
  // InterfaceFlux[0]=InterfaceFlux[0]+max(0,)


#undef LAPIDUS_ARTIFICIAL_VISCOSITY
#ifdef LAPIDUS_ARTIFICIAL_VISCOSITY
  double u1 = oldgrid[ii][jj][0] _MOMX / oldgrid[ii + 1][jj][0] _MASS;
  double u2 = oldgrid[ii - 1][jj][0] _MOMX / oldgrid[ii - 1][jj][0] _MASS;
  double v1 = oldgrid[ii][jj + 1][0] _MOMY / oldgrid[ii][jj + 1][0] _MASS;
  double v2 = oldgrid[ii - 1][jj + 1][0] _MOMY / oldgrid[ii][jj - 1][0] _MASS;
  double v3 = oldgrid[ii][jj - 1][0] _MOMY / oldgrid[ii][jj + 1][0] _MASS;
  double v4 = oldgrid[ii - 1][jj - 1][0] _MOMY / oldgrid[ii][jj - 1][0] _MASS;
  v1 = 0.5 * (v1 + v2);
  v2 = 0.5 * (v3 + v4);
  double divv = 0.5 * (u1 - u2 + v1 - v2);
  double delu = 0;


  for (hh = 0; hh < ne; hh++)
    {
      delu =
	(oldgrid[ii][jj][0].array[hh] -
	 oldgrid[ii - d[0]][jj - d[1]][0].array[hh]);
      InterfaceFlux[hh] =
	InterfaceFlux[hh] + dtodx * 0.1 * min (0.0, divv) * delu;
    }
#endif /* LAPIDUS_ARTIFICIAL_VISCOSITY */

#endif /* ROE */

#ifdef VANLEER
  rc =
    vanleer_flux_vector_split (rightstate, rightstate, fp, fn, ii, jj, idir);
#endif

  delete leftstate;
  delete rightstate;
  delete flux;
  delete fp;
  delete fn;


  return 0;

}

int
ptoc (double *p, double *c)
{
  double rho;
  double px;
  double py, et, velx, vely, velz, velx2, vely2, velz2, ke, pressure;
  double bx, by, bz, bsquared;
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
     cout << "p " << p[0] 
     << " " << p[1] 
     << " " << p[2] 
     << " " << p[3] 
     << endl;
     cout << "c " << c[0] 
     << " " << c[1] 
     << " " << c[2] 
     << " " << c[3] 
     << endl;
   */

  if (c[4] <= 1e-20)
    {
      cout << "p " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " "
	<< p[4] << " " << p[5] << " " << p[6] << " " << p[7] << endl;
      cout << "c " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << " "
	<< c[4] << " " << c[5] << " " << c[6] << " " << c[7] << endl;
      cout << "ptoc: Mein Gott!" << endl;
//              exit (0);
    }

  if (vely < -400)
    {
      cout << setiosflags (ios::fixed);
      cout << "p " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " "
	<< p[4] << " " << p[5] << " " << p[6] << " " << p[7] << endl;
      cout << "c " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << " "
	<< c[4] << " " << c[5] << " " << c[6] << " " << c[7] << endl;
    }

  return 0;
}

int
ctop (double *c, double *p)
{
  double rho, rhoi;
  double px;
  double py, pz, et, velx, vely, velz, velx2, vely2, velz2, ke, pressure;
  double bx, by, bz, bsquared;
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
  bsquared = bx * bx + by * by + bz * bz;
  pressure = et - ke - 0.5 * bsquared;
  pressure = pressure * (gammam1);
  p[0] = rho;
  p[1] = velx;
  p[2] = vely;
  p[3] = velz;
  p[4] = pressure;
  p[5] = bx;
  p[6] = by;
  p[7] = bz;


#ifdef VERBOSE_PRESSURE
  if (pressure < 1e-20)
    {
      cout << "p " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " "
	<< p[4] << " " << p[5] << " " << p[6] << " " << p[7] << endl;
      cout << "c " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << " "
	<< c[4] << " " << c[5] << " " << c[6] << " " << c[7] << endl;
      cout << "ctop: pressure low" << endl;
      if (pressure < 0)
	{
	  cout << "pressure <0 " << pressure << endl;
	  return 1;
	}
    }
#endif
  if (vely < -400)
    {
      cout << setiosflags (ios::fixed);
      cout << "p " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " "
	<< p[4] << " " << p[5] << " " << p[6] << " " << p[7] << endl;
      cout << "c " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << " "
	<< c[4] << " " << c[5] << " " << c[6] << " " << c[7] << endl;
      cout << "Vel < 400" << endl;
    }
  return 0;
}

int
flatten_slopes (Array2D < zone > oldgrid, double a_flat_r, double a_flat_l)
{
  // Flatten slopes of muscl reconstruction before going into the Riemann solver
  double top;
  double bottom;
//  top = fabs ();
  double chi = 0;
  return 0;
}
