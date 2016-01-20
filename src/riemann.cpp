/* $Id: riemann.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "riemann.h"
#include "PhysConsts.h"
#define DEBUG_RIEMANN
#define DEBUG_EV
#undef DEBUG_RIEMANN
//#define CONS_FLUX
#define EINFELDT_FIX
#define MINVAL(a,b) a<b?a:b
#define MAXVAL(a,b) a>b?a:b
/* 
currently no entropy fix is in place 
*/
int
riemann (double *leftstate,
	 double *rightstate,
	 double *roe_flux,
	 double *res_state, int time_step, double *max_speed, int idir)
{
  double gammag = PhysConsts::gamma;

  double pmin = 1e-10;
  ofstream outFile;
  char outfilename[50] = "output/riemann.log";

  int rc = 0;
  int ii = 0;
  int jj = 0;
  int kk = 0;
  int hh = 0;

  double av_state[8];
  double rl = leftstate[0];
  double ul = leftstate[1];
  double vl = leftstate[2];
  double wl = leftstate[3];
  double pl = leftstate[4];
  double bul = leftstate[5];
  double bvl = leftstate[6];
  double bwl = leftstate[7];

  double rr = rightstate[0];
  double ur = rightstate[1];
  double vr = rightstate[2];
  double wr = rightstate[3];
  double pr = rightstate[4];
  double bur = rightstate[5];
  double bvr = rightstate[6];
  double bwr = rightstate[7];
  double cslow = 0;
  double cslow2 = 0;
  double calfven = 0;
  double calfven2 = 0;
  double cfast = 0;
  double cfast2 = 0;
  double bsquared = 0;
  double gammam1 = gammag - 1;
  double gammam1i = 1.0 / gammam1;


  double rho_rl = 0;
  double u_rl = 0;
  double v_rl = 0;
  double w_rl = 0;
  double p_rl = 0;
  double bu_rl = 0;
  double bv_rl = 0;
  double bw_rl = 0;

  double lambda[7];
  double lstate[7];
  double rstate[7];
  double lrsp[7];
  double rrsp[7];
  double fl[7], fr[7], fx[7];
  double eigenwt[7];
  double temp;

  double rho = 0;
  double mass = 0;
  double rhoi = 0;
  double srhoi = 0;
  double u = 0;
  double v = 0;
  double w = 0;
  double bu = 0;
  double bv = 0;
  double bw = 0;
  double p = 0;
  double csound2 = 0;
  double term = 0;
  double a_star2 = 0;
  double vv2 = 0;
  double kinetic = 0;
  double p_magnetic = 0;
  double internal_energy = 0;
  double energy = 0;
  double v2 = 0;
  double vdotb = 0;

  double cc[7];
  double dvy = 0, dvz = 0;

  if (rl < 0 || rr < 0)
    {
      cerr << "Negative density in roe solver!" << endl;
      exit (0);
    }
  outFile.open (outfilename, ofstream::app);
  if (!outFile)
    {
      cerr << "Unable to open file" << endl;
    }

  if (idir == 1)
    {
      /* Keep fluxes */
    }
  else if (idir == 2)
    {
      rl = leftstate[0];
      ul = leftstate[2];
      vl = leftstate[3];
      wl = leftstate[1];
      pl = leftstate[4];
      pl = std::max (pl, pmin);
      bul = leftstate[6];
      bvl = leftstate[7];
      bwl = leftstate[5];

      rr = rightstate[0];
      ur = rightstate[2];
      vr = rightstate[3];
      wr = rightstate[1];
      pr = rightstate[4];
      pr = std::max (pr, pmin);
      bur = rightstate[6];
      bvr = rightstate[7];
      bwr = rightstate[5];
      /* Rotate fluxes */
    }
  else if (idir == 3)
    {
      /* Rotate fluxes */
      rl = leftstate[0];
      ul = leftstate[3];
      vl = leftstate[1];
      wl = leftstate[2];
      pl = leftstate[4];
      bul = leftstate[7];
      bvl = leftstate[5];
      bwl = leftstate[6];

      rr = rightstate[0];
      ur = rightstate[3];
      vr = rightstate[1];
      wr = rightstate[2];
      pr = rightstate[4];
      bur = rightstate[7];
      bvr = rightstate[5];
      bwr = rightstate[6];
    }
  /* Convert conserved to primitives */
  ul = ul / rl;
  vl = vl / rl;
  wl = wl / rl;
  v2 = ul * ul + vl * vl + wl * wl;
  bsquared = bul * bul + bvl * bvl + bwl * bwl;
  pl = pl - (double) 0.5 *(rl * v2) - 0.5 * bsquared;
  pl = pl * gammam1;
  pl = max (pl, pmin);

  ur = ur / rr;
  vr = vr / rr;
  wr = wr / rr;
  v2 = ur * ur + vr * vr + wr * wr;
  bsquared = bur * bur + bvr * bvr + bwr * bwr;
  pr = pr - (double) 0.5 *(rr * v2) - 0.5 * bsquared;
  pr = pr * gammam1;
  pr = max (pr, pmin);
  if (pl < 0)
    {
      cout << __FUNCTION__ << ": Negative pressure in roe solver! p = " << pl
	<< endl;
      cout << " rr " << rr << " ur " << ur << " vr " << vr << " wr " << wr <<
	" pr " << pr << " bur " << bur << " bvr " << bvr << " bwr " << bwr <<
	endl;
      cout << " rl " << rl << " ul " << ul << " vl " << vl << " wl " << wl <<
	" pl " << pl << " bul " << bul << " bvl " << bvl << " bwl " << bwl <<
	endl;
      return (1);
    }
  if (pr < 0)
    {
      cout << __FUNCTION__ << ": Negative pressure in roe solver! p = " << pr
	<< endl;
      cout << " rr " << rr << " ur " << ur << " vr " << vr << " wr " << wr <<
	" pr " << pr << " bur " << bur << " bvr " << bvr << " bwr " << bwr <<
	endl;
      cout << " rl " << rl << " ul " << ul << " vl " << vl << " wl " << wl <<
	" pl " << pl << " bul " << bul << " bvl " << bvl << " bwl " << bwl <<
	endl;
      return (1);
    }

  lrsp[0] = rl;
  lrsp[1] = ul;
  lrsp[2] = vl;
  lrsp[3] = wl;
  lrsp[4] = pl;
  lrsp[5] = bvl;
  lrsp[6] = bwl;

  rrsp[0] = rr;
  rrsp[1] = ur;
  rrsp[2] = vr;
  rrsp[3] = wr;
  rrsp[4] = pr;
  rrsp[5] = bvr;
  rrsp[6] = bwr;

  for (ii = 0; ii < 7; ii++)
    {
      lstate[ii] = lrsp[ii];
      rstate[ii] = rrsp[ii];
    }
  cc[0] = (rr - rl);
  cc[1] = (ur - ul);
  cc[2] = (vr - vl);
  dvy = vr - vl;
  cc[3] = (wr - wl);
  dvz = wr - wl;
  cc[4] = (pr - pl);
  if (fabs (cc[4]) < 1e-10)
    cc[4] = 0;
  cc[5] = (bvr - bvl);
  cc[6] = (bwr - bwl);

  /* compute the averaged quantities */
  rho_rl = 0.5 * (rr + rl);
  u_rl = 0.5 * (ul + ur);
  if (u_rl == 0)
    {
      /* work out res state and flux and exit? */
    }
  v_rl = 0.5 * (vl + vr);
  w_rl = 0.5 * (wl + wr);
  p_rl = 0.5 * (pl + pr);
  bu_rl = 0.5 * (bul + bur);
  bv_rl = 0.5 * (bvl + bvr);
  bw_rl = 0.5 * (bwl + bwr);

  av_state[0] = rho_rl;
  av_state[1] = u_rl;
  av_state[2] = v_rl;
  av_state[3] = w_rl;
  av_state[4] = p_rl;
  av_state[5] = bu_rl;
  av_state[6] = bv_rl;
  av_state[7] = bw_rl;

  /* Allocte memory for eigenvectors */
  Array2D < double >levec (7, 7);
  Array2D < double >revec (7, 7);
  Array2D < double >levc (7, 7);
  Array2D < double >revc (7, 7);

  /* Compute fast and slow speeds */
  rho = rho_rl;
  rhoi = 1 / rho;
  srhoi = 1 / sqrt (rho);
  u = u_rl;
  v = v_rl;
  w = w_rl;
  bu = bu_rl * sqrt (rhoi);
  bv = bv_rl * sqrt (rhoi);
  bw = bw_rl * sqrt (rhoi);
  bsquared = (bu * bu + bv * bv + bw * bw);
  vv2 = (u * u + v * v + w * w);
  p = p_rl;
  calfven2 = fabs (bu * bu);
  calfven = sqrt (calfven2);
  csound2 = gammag * p * rhoi;
  a_star2 = csound2 + bsquared;
  term = sqrt (a_star2 * a_star2 - 4.0 * csound2 * calfven2);
  term = fabs (term);
  cslow2 = 0.5 * (a_star2 - term);
  if (cslow2 < 0)
    {
      if (cslow2 < -1e-10)
	{
	  cout << "Slow speed went negative!! " << cslow2 << endl;
	  cout
	    << " csound2 " << csound2
	    << " calfven2 " << calfven2
	    << " bsquared " << bsquared
	    << " bx " << bu << " by " << bv << " bz " << bw << endl;
	  return 1;
	}
      else
	{
	  cslow2 = 0;
	}
    }
  cfast2 = 0.5 * (a_star2 + term);
  cslow = sqrt (cslow2);
  cfast = sqrt (cfast2);

  /* compute left and right characteristic eigenvectors */
  rc = eigenvectors (av_state, levec, revec, levc, revc, dvy, dvz);

  /* compute eigenvalues */
  lambda[0] = u_rl - cfast;
  lambda[1] = u_rl - calfven;
  lambda[2] = u_rl - cslow;
  lambda[3] = u_rl;
  lambda[4] = u_rl + cslow;
  lambda[5] = u_rl + calfven;
  lambda[6] = u_rl + cfast;


#ifdef DEBUG_RIEMANN
  if (ul != ur)
    {
      outFile << endl;
      outFile << "Time Step " << time_step << endl;
      outFile << "Left and Right States -----" << endl;

      outFile << "rl " << rl
	<< " ul " << ul
	<< " vl " << vl
	<< " wl " << wl
	<< " pl " << pl
	<< " bul " << bul << " bvl " << bvl << " bwl " << bwl << endl;

//               for ( hh=0 ; hh<8 ; hh++ )
//               {
//                       outFile << lrsp[hh] << " " ;
//               }
//               outFile << endl;



      outFile << "rr " << rr
	<< " ur " << ur
	<< " vr " << vr
	<< " wr " << wr
	<< " pr " << pr
	<< " bur " << bur << " bvr " << bvr << " bwr " << bwr << endl;
      outFile << endl;

      outFile << "Average State -----" << endl;
      outFile << "rl " << rho_rl
	<< " ul " << u_rl
	<< " vl " << v_rl
	<< " wl " << w_rl
	<< " pl " << p_rl
	<< " bul " << bu_rl << " bvl " << bv_rl << " bwl " << bw_rl << endl;

      outFile << endl;
      outFile << "Difference Between Left and Right States -----" << endl;
      outFile << "del_rho " << cc[0]
	<< " del_u " << cc[1]
	<< " del_v " << cc[2]
	<< " del_w " << cc[3]
	<< " del_p " << cc[4]
	<< " del_bv " << cc[5] << " del_bw " << cc[6] << endl;
      outFile << endl;

      outFile << "Fast, Slow , Alfven speeds-----" << endl;
      outFile
	<< " cfast " << cfast
	<< " cslow " << cslow << " calfven " << calfven << endl;
      outFile << endl;

      outFile << "Characteristic wave speeds (Eigenvalues)-----" << endl;
      outFile << "lam1 " << lambda[0]
	<< " lam2 " << lambda[1]
	<< " lam3 " << lambda[2]
	<< " lam4 " << lambda[3]
	<< " lam5 " << lambda[4]
	<< " lam6 " << lambda[5] << " lam7 " << lambda[6] << endl;
      outFile << endl;

#ifdef DEBUG_EV



      outFile << "Left Eigenvectors -----" << endl;
      for (jj = 0; jj < 7; jj++)
	{
	  for (ii = 0; ii < 7; ii++)
	    {
	      outFile << levec[ii][jj] << "\t";
	    }
	  outFile << endl;
	}
      outFile << endl;


      outFile << "Right Eigenvectors -----" << endl;
      for (jj = 0; jj < 7; jj++)
	{
	  for (ii = 0; ii < 7; ii++)
	    {
	      outFile << revec[ii][jj] << "\t";
	    }
	  outFile << endl;
	}
      outFile << endl;



      outFile << "Left Cons Eigenvectors -----" << endl;
      for (jj = 0; jj < 7; jj++)
	{
	  for (ii = 0; ii < 7; ii++)
	    {
	      outFile << levc[ii][jj] << "\t";
	    }
	  outFile << endl;
	}
      outFile << endl;


      outFile << "Right Cons Eigenvectors -----" << endl;
      for (jj = 0; jj < 7; jj++)
	{
	  for (ii = 0; ii < 7; ii++)
	    {
	      outFile << revc[ii][jj] << "\t";
	    }
	  outFile << endl;
	}
      outFile << endl;

#endif /* DEBUG_EV */

      outFile << "Dots of left and right ev -----" << endl;
      for (kk = 0; kk < 7; kk++)
	{
	  temp = 0;
	}
      for (kk = 0; kk < 7; kk++)
	{
	  for (jj = 0; jj < 7; jj++)
	    {
	      outFile << levec[jj][kk] * revec[jj][kk] << " ";
	      temp = temp + levec[jj][kk] * revec[jj][kk];
	    }
	  outFile << temp << endl;
	  temp = 0;
	}


      outFile << "Dots of left and right ev -----" << endl;
      for (kk = 0; kk < 7; kk++)
	{
	  temp = 0;
	}
      for (kk = 0; kk < 7; kk++)
	{
	  for (jj = 0; jj < 7; jj++)
	    {
	      outFile << levc[jj][kk] * revc[jj][kk] << " ";
	      temp = temp + levc[jj][kk] * revc[jj][kk];
	    }
	  outFile << temp << endl;
	  temp = 0;
	}
    }



#endif /* DEBUG_RIEMANN */

  /* Compute resolved state primitives */
#ifdef DEBUG_RIEMANN
  if (ul != ur)
    {
      outFile << endl;
    }
  if (ul != ur)
    {
      outFile << "Eigenweights-----" << endl;
    }
#endif /* DEBUG_RIEMANN */
  for (jj = 0; jj < 7; jj++)
    {
      eigenwt[jj] = 0;
      for (kk = 0; kk < 7; kk++)
	{
	  eigenwt[jj] = eigenwt[jj] + levec[kk][jj] * cc[kk];
#ifdef DEBUG_RIEMANN
	  if (ul != ur)
	    {
	      outFile << eigenwt[jj] << "\t\t";
	    }
#endif /* DEBUG_RIEMANN */
	}
#ifdef DEBUG_RIEMANN
      if (ul != ur)
	{
	  outFile << endl;
	}
#endif /* DEBUG_RIEMANN */
    }
#ifdef DEBUG_RIEMANN
  if (ul != ur)
    {
      outFile << endl;
    }
  if (ul != ur)
    {
      outFile << "Left Primitive State-----" << endl;
    }
#endif /* DEBUG_RIEMANN */


  for (jj = 0; jj < 7; jj++)
    {
      for (kk = 0; kk < 7; kk++)
	{
	  if (lambda[kk] < 0)
	    {
	      lrsp[jj] = lrsp[jj] + eigenwt[kk] * revec[jj][kk];
	    }
#ifdef DEBUG_RIEMANN
	  if (ul != ur)
	    {
	      outFile << lrsp[jj] << "\t\t";
	    }
#endif /* DEBUG_RIEMANN */
	}
#ifdef DEBUG_RIEMANN
      if (ul != ur)
	{
	  outFile << endl;
	}
#endif /* DEBUG_RIEMANN */
    }
  bu = bul;


#ifdef DEBUG_RIEMANN
  if (ul != ur)
    {
      outFile << endl;
    }
  if (ul != ur)
    {
      outFile << "Right Primitive State-----" << endl;
    }
#endif /* DEBUG_RIEMANN */

  for (jj = 0; jj < 7; jj++)
    {
//                        rrsp[jj]=rrsp[jj];
      for (kk = 6; kk > -1; kk--)
	{
	  if (lambda[kk] > 0)
	    {
	      rrsp[jj] = rrsp[jj] - eigenwt[kk] * revec[jj][kk];

	    }
#ifdef DEBUG_RIEMANN
	  if (ul != ur)
	    {
	      outFile << rrsp[jj] << "\t\t";
	    }
#endif /* DEBUG_RIEMANN */
	  bu = bur;
	}

#ifdef DEBUG_RIEMANN
      if (ul != ur)
	{
	  outFile << endl;
	}
#endif /* DEBUG_RIEMANN */
    }


  if (lambda[3] > 1e-10)
    {
      mass = lrsp[0];
      u = lrsp[1];
      v = lrsp[2];
      w = lrsp[3];
      p = lrsp[4];
      bu = bul;
      bv = lrsp[5];
      bw = lrsp[6];
    }
  else if (lambda[3] < -1e-10)
    {
      mass = rrsp[0];
      u = rrsp[1];
      v = rrsp[2];
      w = rrsp[3];
      p = rrsp[4];
      bu = bur;
      bv = rrsp[5];
      bw = rrsp[6];
    }
  else
    {
      mass = 0.5 * (lrsp[0] + rrsp[0]);
      u = 0.5 * (lrsp[1] + rrsp[1]);
      v = 0.5 * (lrsp[2] + rrsp[2]);
      w = 0.5 * (lrsp[3] + rrsp[3]);
      p = 0.5 * (lrsp[4] + rrsp[4]);
      bu = bur;
      bv = 0.5 * (lrsp[5] + rrsp[5]);
      bw = 0.5 * (lrsp[6] + rrsp[6]);

      /*
         mass = 0.5*(lstate[0]+rstate[0]); 
         u    = 0.5*(lstate[1]+rstate[1]);
         v    = 0.5*(lstate[2]+rstate[2]);
         w    = 0.5*(lstate[3]+rstate[3]);
         p    = 0.5*(lstate[4]+rstate[4]);
         bu   = bur;
         bv   = 0.5*(lstate[5]+rstate[5]);
         bw   = 0.5*(lstate[6]+rstate[6]);
       */
    }


#ifdef DEBUG_RIEMANN
  if (ul != ur)
    {
      outFile << endl;
      outFile << "Resolved state primitives -----" << endl;
      outFile << "mass " << mass << " u " << u <<
	" v " << v <<
	" w " << w <<
	" p " << p << " bu " << bu << " bv " << bv << " bw " << bw << endl;
      outFile << endl;
    }
#endif /* DEBUG_RIEMANN */

#ifdef CONS_FLUX

  mass = lstate[0];
  u = lstate[1];
  v = lstate[2];
  w = lstate[3];
  p = lstate[4];
  bu = bul;
  bv = lstate[5];
  bw = lstate[6];


  v2 = u * u + v * v + w * w;
  kinetic = 0.5 * mass * v2;
  bsquared = bu * bu + bv * bv + bw * bw;
  p_magnetic = 0.5 * bsquared;
  internal_energy = p * gammam1i;
  energy = kinetic + p_magnetic + internal_energy;
  vdotb = u * bu + v * bv + w * bw;


  fl[0] = mass * u;
  fl[1] = mass * u * u + p + p_magnetic - bu * bu;
  fl[2] = mass * u * v - bu * bv;
  fl[3] = mass * u * w - bu * bw;
  fl[4] = (energy + p + p_magnetic) * u - bu * (vdotb);
  fl[5] = u * bv - v * bu;
  fl[6] = u * bw - w * bu;


  mass = rstate[0];
  u = rstate[1];
  v = rstate[2];
  w = rstate[3];
  p = rstate[4];
  bu = bul;
  bv = rstate[5];
  bw = rstate[6];


  v2 = u * u + v * v + w * w;
  kinetic = 0.5 * mass * v2;
  bsquared = bu * bu + bv * bv + bw * bw;
  p_magnetic = 0.5 * bsquared;
  internal_energy = p * gammam1i;
  energy = kinetic + p_magnetic + internal_energy;
  vdotb = u * bu + v * bv + w * bw;


  fr[0] = mass * u;
  fr[1] = mass * u * u + p + p_magnetic - bu * bu;
  fr[2] = mass * u * v - bu * bv;
  fr[3] = mass * u * w - bu * bw;
  fr[4] = (energy + p + p_magnetic) * u - bu * (vdotb);
  fr[5] = u * bv - v * bu;
  fr[6] = u * bw - w * bu;

#ifdef DEBUG_RIEMANN
  if (ul != ur)
    {
      for (ii = 0; ii < 7; ii++)
	outFile << " " << fl[ii];
      outFile << endl;
      for (ii = 0; ii < 7; ii++)
	outFile << " " << fr[ii];
      outFile << endl;
    }
#endif

#ifdef EINFELDT_FIX
  {
    double eval_l[7];
    double eval_r[7];
    rc = eigenvalues (lstate, eval_l);
    rc = eigenvalues (rstate, eval_r);
    double bplus = MAXVAL (eval_r[6], lambda[6]);
    double bminus = MINVAL (eval_l[6], lambda[6]);
    bplus = MAXVAL (bplus, 0.0);
    bminus = MINVAL (bminus, 0.0);
    double denom = cfast + 0.5 * fabs (MAXVAL (eval_r[6], lambda[6]) +
				       MINVAL (eval_l[6], lambda[6]));
    double delta = cfast / denom;
    double coef0 = 1 / (bplus - bminus);
    double coef1 = (bplus + bminus) * coef0;
    double coef2 = (bplus * bminus) * coef0 * 2.0;
    lambda[0] = coef1 * lambda[0] - coef2;
    lambda[1] = coef1 * lambda[1] - coef2 * (1.0 - delta);
    lambda[2] = coef1 * lambda[2] - coef2;
    lambda[3] = coef1 * lambda[3] - coef2 * (1.0 - delta);
    lambda[4] = coef1 * lambda[4] - coef2;
    lambda[5] = coef1 * lambda[5] - coef2 * (1.0 - delta);
    lambda[6] = coef1 * lambda[6] - coef2;
  }
// entropy fix 
//      lambda[kk]=0;
#endif /* EINFELDT_FIX */


  for (ii = 0; ii < 7; ii++)
    {
      fx[ii] = 0.5 * (fl[ii] + fr[ii]);
#ifdef DEBUG_RIEMANN
      if (ul != ur)
	{
	  outFile << " " << fx[ii];
	}
#endif
    }
#ifdef DEBUG_RIEMANN
  if (ul != ur)
    {
      outFile << endl;
    }
#endif
  for (jj = 0; jj < 7; jj++)
    {
      for (kk = 0; kk < 7; kk++)
	{
	  fx[jj] =
	    fx[jj] - 0.5 * fabs (lambda[kk]) * eigenwt[kk] * revc[jj][kk];
#ifdef DEBUG_RIEMANN
	  if (ul != ur)
	    {
	      outFile << " " << fx[jj];
	    }
#endif
	}
#ifdef DEBUG_RIEMANN
      if (ul != ur)
	{
	  outFile << endl;
	}
#endif

    }

#endif

  v2 = u * u + v * v + w * w;
  kinetic = 0.5 * mass * v2;
  bsquared = bu * bu + bv * bv + bw * bw;
  p_magnetic = 0.5 * bsquared;
  internal_energy = p * gammam1i;
  energy = kinetic + p_magnetic + internal_energy;
  vdotb = u * bu + v * bv + w * bw;



  /* compute the fluxes */
  roe_flux[0] = mass * u;
  roe_flux[1] = mass * u * u + p + p_magnetic - bu * bu;
  roe_flux[2] = mass * u * v - bu * bv;
  roe_flux[3] = mass * u * w - bu * bw;
  roe_flux[4] = (energy + p + p_magnetic) * u - bu * (vdotb);
  roe_flux[5] = 0;
  roe_flux[6] = u * bv - v * bu;
  roe_flux[7] = u * bw - w * bu;
  /* Compute resolved state flux from P* */

//      if (delta_rho >= 0.0000001)


  if (idir == 1)
    {
      /* Do nothing */
      res_state[0] = mass;
      res_state[1] = mass * u;
      res_state[2] = mass * v;
      res_state[3] = mass * w;
      res_state[4] = energy;
      res_state[5] = bu;
      res_state[6] = bv;
      res_state[7] = bw;
    }
  else if (idir == 2)
    {
      /* Rotate the fluxes */
      roe_flux[0] = mass * u;
      roe_flux[2] = mass * u * u + p + p_magnetic - bu * bu;
      roe_flux[3] = mass * u * v - bu * bv;
      roe_flux[1] = mass * u * w - bu * bw;
      roe_flux[4] = (energy + p + p_magnetic) * u - bu * (vdotb);
      roe_flux[6] = 0;
      roe_flux[7] = u * bv - v * bu;
      roe_flux[5] = u * bw - w * bu;

      res_state[0] = mass;
      res_state[2] = mass * u;
      res_state[3] = mass * v;
      res_state[1] = mass * w;
      res_state[4] = energy;
      res_state[6] = bu;
      res_state[7] = bv;
      res_state[5] = bw;
    }
  else if (idir == 3)
    {
      /* Rotate the fluxes */
      roe_flux[0] = mass * u;
      roe_flux[1] = mass * u * w - bu * bw;
      roe_flux[2] = mass * u * v - bu * bv;
      roe_flux[3] = mass * u * u + p + p_magnetic - bu * bu;
      roe_flux[4] = (energy + p + p_magnetic) * u - bu * (vdotb);
      roe_flux[5] = u * bv - v * bu;
      roe_flux[6] = u * bw - w * bu;
      roe_flux[7] = 0;
    }


#ifdef DEBUG_RIEMANN
  if (ul != ur)
    {
      outFile << "Flux Vector -----" << endl;
      outFile << roe_flux[0]
	<< " " << roe_flux[1]
	<< " " << roe_flux[2]
	<< " " << roe_flux[3]
	<< " " << roe_flux[4]
	<< " " << roe_flux[5]
	<< " " << roe_flux[6] << " " << roe_flux[7] << endl;
      outFile << "##############" << endl;
      outFile << endl;
    }
#endif /* DEBUG_RIEMANN */

#ifdef CONS_FLUX
  if (idir == 1)
    {
      for (kk = 0; kk < 7; kk++)
	{
	  roe_flux[kk] = fx[kk];
	}
      roe_flux[5] = 0;
      roe_flux[6] = fx[5];
      roe_flux[7] = fx[6];
    }
  else if (idir == 2)
    {

      roe_flux[0] = fx[0];
      roe_flux[1] = fx[2];
      roe_flux[2] = fx[3];
      roe_flux[3] = fx[1];
      roe_flux[4] = fx[4];
      roe_flux[5] = fx[6];
      roe_flux[6] = fx[7];
      roe_flux[7] = fx[5];
    }
#endif
  for (ii = 0; ii < 8; ii++)
    {
      if (fabs (roe_flux[ii]) < -1e10)
	{
	  roe_flux[ii] = 0;
	}
    }

  outFile.close ();
  return 0;
}
