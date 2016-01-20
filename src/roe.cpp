/* $Id: roe.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "roe.h"
#define SERIOUS_LOGGING
#undef SERIOUS_LOGGING
int
roe (double *leftstate, double *rightstate, double *roeflux, double *maxspeed,
     int iii, int jjj, int timestep, int idir)
{


  ofstream outFile;
  char outfilename[50] = "output/roe.log";
  int ii, jj;
  int d[2];

  double gammai = 1.0 / gammag;
  double gammam1 = gammag - 1;
  double gammam1i = 1.0 / gammam1;
  double rr, ri, px, py, et, ke;
  double rl, ul, vl, pl, hl, al;
  double ur, vr, pr, hr, ar;
  double rho_rl, url, vrl, prl, arl, hrl;
  double lambda[ne];
  double kx = 0, ky = 0;
  double eigenwt[ne];
  double rev[ne][ne];
  double cc[ne];
  double lres_state_prim[ne];
  double rres_state_prim[ne];
  double smallp = 1e-5;
  double rho_i_2_c = 0;
  double delta_w[ne];
  double delta_p;
  double delta_v;
  double delta_u;
  double delta_rho;
  double rflux[ne];
  double lflux[ne];

//  cout << " Roe solve" << endl;
  if (idir == 1)
    {
      kx = 1.0;
      ky = 0.0;
      d[0] = 1;
      d[1] = 0;
    }
  else if (idir == 2)
    {
      kx = 0.0;
      ky = 1.0;
      d[0] = 0;
      d[1] = 1;
    }
  else
    {
      exit (0);
    }

  /* Determine the x fluxes */
  /* Linear interpolate to get left and right states */
  /* Conserved left state */

  rl = leftstate[0];
  px = leftstate[1];
  py = leftstate[2];
  et = leftstate[3];

  ri = 1.0 / rl;
  ul = px * ri;
  vl = py * ri;
  ke = 0.5 * rl * (ul * ul + vl * vl);
  pl = et - ke;
  pl = pl * gammam1;
#ifdef DEBUG
  if (pl < 0.0)
    {
      cout << "roe: Negative pressure at " << iii << "," << jjj << endl;
      cout
	<< "rl= " << rl
	<< " ul= " << ul
	<< " vl= " << vl << " et= " << et << " pl= " << pl << endl;
      exit (0);
      pl = smallp;
    }
#endif /* DEBUG */

  lflux[0] = rl * ul;
  lflux[1] = rl * ul * ul + pl;
  lflux[2] = rl * ul * vl;
  lflux[3] = ul * (0.5 * rl * (ul * ul + ul * vl) + pl + pl * gammam1i);

  /* Initialising left resolved state primitives with the left state */
  lres_state_prim[0] = rl;
  lres_state_prim[1] = ul;
  lres_state_prim[2] = vl;
  lres_state_prim[3] = pl;
  hl = ri * (et + pl);
  al = sqrt (gammag * pl / rl);
  /* Conserved right state */


  rr = rightstate[0];
  px = rightstate[1];
  py = rightstate[2];
  et = rightstate[3];



  ri = 1.0 / rr;
  ur = px * ri;
  vr = py * ri;
  ke = 0.5 * rr * (ur * ur + vr * vr);
  pr = et - ke;
  pr = pr * gammam1;
  /* Initialising right resolved state primitives with the right state */
  rres_state_prim[0] = rr;
  rres_state_prim[1] = ur;
  rres_state_prim[2] = vr;
  rres_state_prim[3] = pr;
  if (pr < 0.0)
    {
      cout << "roe: Negative pressure at " << iii << "," << jjj << endl;
      cout
	<< " Rho= " << rr
	<< " Momx= " << px
	<< " Momy= " << py
	<< " Velx= " << ur
	<< " Vely= " << vr << " Etot= " << et << " Pres= " << pr << endl;
      exit (0);
      pr = smallp;
    }
  hr = ri * (et + pr);
  ar = sqrt (gammag * pr / rr);


  rflux[0] = rr * ur;
  rflux[1] = rr * ur * vr;
  rflux[2] = rr * ur * ur + pr;
  rflux[3] = ur * (0.5 * rr * (ur * ur + vr * vr) + pr + pr * gammam1i);

  cc[0] = rr - rl;
  cc[1] = ur - ul;
  cc[2] = vr - vl;
  cc[3] = pr - pl;




  delta_rho = (rr - rl);
  delta_p = (pr - pl);
  delta_u = (ur - ul);
  delta_v = (vr - vl);

  /* Arithmetic mean */
  rho_rl = 0.5 * (rl + rr);
  url = 0.5 * (ul + ur);
  vrl = 0.5 * (vl + vr);
  prl = 0.5 * (pl + pr);

  /* Compute Roe-average state */
  rho_rl = sqrt (rl * rr);
  url = (sqrt (rl) * ul + sqrt (rr) * ur) / (sqrt (rl) + sqrt (rr));
  vrl = (sqrt (rl) * vl + sqrt (rr) * vr) / (sqrt (rl) + sqrt (rr));
  hrl = (sqrt (rl) * hl + sqrt (rr) * hr) / (sqrt (rl) + sqrt (rr));
  arl = sqrt (gammam1 * (hrl - 0.5 * (url * url + vrl * vrl)));
  prl = rho_rl * arl * arl * gammai;




  arl = sqrt (gammag * prl / rho_rl);

  /* Compute wave speeds / eigenvalues */
  if (idir == 1)
    {



      lflux[0] = rl * ul;
      lflux[1] = rl * ul * ul + pl;
      lflux[2] = rl * ul * vl;
      lflux[3] = ul * (0.5 * rl * (ul * ul + vl * vl) + pl + pl * gammam1i);

      lambda[0] = url;
      lambda[1] = url;
      lambda[2] = url + arl;
      lambda[3] = url - arl;

      rflux[0] = rr * ur;
      rflux[1] = rr * ur * ur + pr;
      rflux[2] = rr * ur * vr;
      rflux[3] = ur * (0.5 * rr * (ur * ur + vr * vr) + pr + pr * gammam1i);

    }
  else if (idir == 2)
    {

      lflux[0] = rl * vl;
      lflux[1] = rl * ul * vl;
      lflux[2] = rl * vl * vl + pl;
      lflux[3] = vl * (0.5 * rl * (ul * ul + vl * vl) + pl + pl * gammam1i);

      rflux[0] = rr * vr;
      rflux[1] = rr * ur * vr;
      rflux[2] = rr * vr * vr + pr;
      rflux[3] = vr * (0.5 * rr * (ur * ur + vr * vr) + pr + pr * gammam1i);


      lambda[0] = vrl;
      lambda[1] = vrl;
      lambda[2] = vrl + arl;
      lambda[3] = vrl - arl;
    }

  /* compute the wave strengths */

  delta_w[0] = delta_rho - (delta_p / pow (arl, 2));
  delta_w[1] = ky * delta_u - kx * delta_v;
  delta_w[2] = kx * delta_u + ky * delta_v + (delta_p / (rho_rl * arl));
  delta_w[3] = -kx * delta_u - ky * delta_v + (delta_p / (rho_rl * arl));


  /* Compute Eigenvectors */


  rev[0][0] = 1.0;
  rev[0][1] = url;
  rev[0][2] = vrl;
  rev[0][3] = 0.5 * (url * url + vrl * vrl);

  rev[1][0] = 0.0;
  rev[1][1] = rho_rl * ky;
  rev[1][2] = rho_rl * (-kx);
  rev[1][3] = rho_rl * (url * ky - vrl * kx);

  rho_i_2_c = 0.5 * rho_rl / arl;

  rev[2][0] = rho_i_2_c;
  rev[2][1] = rho_i_2_c * (url + arl * kx);
  rev[2][2] = rho_i_2_c * (vrl + arl * ky);
  rev[2][3] = rho_i_2_c * (hrl + arl * (url * kx + vrl * ky));

  rho_i_2_c = 0.5 * rho_rl / arl;
  rev[3][0] = (rho_i_2_c);
  rev[3][1] = (rho_i_2_c) * (url - arl * kx);
  rev[3][2] = (rho_i_2_c) * (vrl - arl * ky);
  rev[3][3] = (rho_i_2_c) * (hrl - arl * (url * kx + vrl * ky));

#ifdef DEBUG_EV

  cout << "Right ev " << endl;
  for (ii = 0; ii < ne; ii++)
    {
      for (jj = 0; jj < ne; jj++)
	{
	  cout << rev[jj][ii] << " ";
	}
      cout << endl;
    }

#endif
  /* Compute eigenweights */
  for (ii = 0; ii < ne; ii++)
    {
      eigenwt[ii] = 0;
      for (jj = 0; jj < ne; jj++)
	{
	  eigenwt[ii] =
	    eigenwt[ii] + fabs (lambda[jj]) * rev[jj][ii] * delta_w[jj];
	}
    }


  for (jj = 0; jj < ne; jj++)
    {
      roeflux[jj] = 0.5 * (lflux[jj] + rflux[jj]) - 0.5 * eigenwt[jj];
      if (isnan (roeflux[jj]))
	{

	  cout << "roeflux[" << jj << "] go crazy! " << endl;
	  cout << " rl " << rl
	    << " ul " << ul
	    << " vl " << vl << " pl " << pl << " al " << al << endl;
	  cout << " rr " << rr
	    << " ur " << ur
	    << " vr " << vr << " pr " << pr << " ar " << ar << endl;
	}
#ifdef DEBUG_fLUX
      if (rl != rr)
	{
	  cout << "roeflux  "
	    << " " << roeflux[jj]
	    << " " << lflux[jj]
	    << " " << rflux[jj] << " " << eigenwt[jj] << endl;
	}
#endif /* DEBUG */
    }



#ifdef SERIOUS_LOGGING

  //                      if ( jjj ==5) {
//                              if ( iii ==5) {
  if (rl != rr)
    {
      outFile.open (outfilename, ofstream::out | ofstream::app);
      if (!outFile)
	{

	  cerr << "Unable to open file" << endl;
	}
      outFile << "Timestep= " << timestep << endl;
      outFile << endl;
      outFile << " ii= " << iii
	<< " jj= " << jjj << " idir= " << idir << endl;
      outFile << "Left state variables " << endl;
      outFile
	<< " rl " << rl
	<< " ul " << ul
	<< " vl " << vl << " pl " << pl << " al " << al << endl;

      outFile << "Right state variables " << endl;
      outFile
	<< " rr " << rr
	<< " ur " << ur
	<< " vr " << vr << " pr " << pr << " ar " << ar << endl;
      outFile << endl;

      outFile
	<< "rr " << rr
	<< " hr " << hr << " rl " << rl << " hl " << hl << endl;

      outFile << "del_rho " << delta_rho
	<< " del_p " << delta_p
	<< " del_u " << delta_u << " del_v " << delta_v << endl;

      outFile << "Roe Averages" << endl;
      outFile << "rho_rl " << rho_rl
	<< " url " << url
	<< " vrl " << vrl << " hrl " << hrl << " arl " << arl << endl;
      outFile << endl;

      outFile << "Characteristic wave speeds " << endl;
      outFile << "lam1 " << lambda[0]
	<< " lam2 " << lambda[1]
	<< " lam3 " << lambda[2] << " lam4 " << lambda[3] << endl;
      outFile << endl;

      outFile << "Characteristic wave strengths " << endl;
      outFile << "del_v1 " << delta_w[0]
	<< " del_v2 " << delta_w[1]
	<< " del_v3 " << delta_w[2] << " del_v4 " << delta_w[3] << endl;
      outFile << endl;

      outFile << "Right Eigenvectors " << endl;
      for (ii = 0; ii < ne; ii++)
	{
	  for (jj = 0; jj < ne; jj++)
	    {
	      outFile << rev[jj][ii] << " ";
	    }
	  outFile << endl;
	}
      outFile << endl;

      for (ii = 0; ii < ne; ii++)
	{
	  for (jj = 0; jj < ne; jj++)
	    {
	      outFile << fabs (lambda[jj]) * delta_w[jj] * rev[jj][ii] << " ";
	    }
	  outFile << endl;
	}

      outFile << endl;
      outFile << "#### Roe Fluxes ####" << endl;
      for (jj = 0; jj < ne; jj++)
	{
	  outFile
	    << "roe_flux [" << jj << "] "
	    << roeflux[jj]
	    << " " << 0.5 * (lflux[jj] + rflux[jj])
	    << " " << 0.5 * eigenwt[jj] << endl;
	}

      outFile << endl;
      outFile.close ();
    }
//                              }

#endif /* SERIOUS_LOGGING */

  return 0;

}
