/* $Id: molcool.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "molcool.h"


int
molcool (double temperature, double n_H2, double nh, double *elossrate)
{

  double Lvh = 0;
  double Kh = 0;
  double Lvl = 0;
  double Lrh = 0;
  double Lrl = 0;
  double xxx = 0;
  double qqq = 0;

#define DEBUG1
#undef DEBUG1
#ifdef DEBUG1
  cout << "molcool: "
    << " " << temperature << " " << n_H2 << " " << nh << endl;
#endif /* DEBUG */

  if (temperature > 1635.0)
    {
      Kh = (1.0e-12) * sqrt (temperature) * exp (-1000.0 / temperature);
    }
  else
    {
      Kh =
	1.4e-13 * exp ((temperature / 125.0) -
		       (temperature / 577.0) * (temperature / 577.0));
    }

  Lvh = (1.1e-13) * exp (-6744.0 / temperature);
  Lvl = (8.18e-13) * (nh * Kh);
  xxx = log10 (temperature / 10000.0);

  if (temperature > 1087.0)
    {
      Lrh = (3.9e-19) * exp (-6118.0 / temperature);
    }
  else
    {
      Lrh = pow (10.0, -19.24 + 0.474 * xxx - 1.274 * xxx * xxx);
    }

  qqq = pow (n_H2, 0.77) + 1.2 * pow (nh, 0.77);
  if (temperature > 4031.0)
    {
      Lrl = (1.38e-22) * qqq * exp (-9243.0 / temperature);
    }
  else
    {
      Lrl = pow (10.0, -22.9 - 0.553 * xxx - 1.148 * xxx * xxx) * qqq;
    }

  xxx = n_H2 * (Lvh / (1 + (Lvh / Lvl)) + Lrh / (1 + (Lrh / Lrl)));
  *elossrate = n_H2 * (Lvh / (1 + (Lvh / Lvl)) + Lrh / (1 + (Lrh / Lrl)));

  if (xxx > 1.0e+06)
    {
      cout << "The h2 cooling is too large" << endl;;
//      exit (0);
    }

  if (xxx < 0.0)
    {
      cout << "The h2 cooling is negative" << endl;;
//      exit (0);
    }

#ifdef DEBUG1
  cout << "exit" << endl;
#endif /* DEBUG1 */
  return 0;
}
