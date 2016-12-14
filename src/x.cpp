/* $Id: x.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
// Energy Term (Joule Heating Energy)
//
//  -B x \eta J
//   \eta is a diagonal tensor or rank 2


int test() {
  double bR = 0;
  double bZ = 0;
  double bPhi = 0;

  //  etam = Va *H etc...
  double etam = alpham * Va * H * exp (-2 * zz * zz / (H * H));
  double chim = 1.0;
  double etamdash = chim * etam;
  double etaR = etamdash;
  double etaZ = etamdash;
  double etaPhi = etam;

  double jR = 0;
  double jZ = 0;
  double jPhi = 0;

  double b2 = 0;
  double b1 = 0;
  // xflux current
  jR =
    0.5 * idy * ((Bz[i][j + 1] + Bz[i - 1][j + 1]) -
		 (Bz[i][j - 1] + Bz[i - 1][j - 1]));
  jZ = (-idx) * ((Bz[i][j] - Bz[i - 1][j]));
  jPhi =
    (idx) * ((By[i][j] - By[i - 1][j]))
    -
    0.5 * idy * ((Bx[i][j + 1] + Bx[i - 1][j + 1]) -
		 (Bx[i][j - 1] + Bx[i - 1][j - 1]));

  // yflux current
  jR = idy * (Bz[i][j] - Bz[i][j - 1]);
  jZ =
    0.5 * (-idx) * ((Bz[i + 1][j] + Bz[i + 1][j - 1]) -
		    (Bz[i - 1][j] + Bz[i - 1][j - 1]));
  jPhi =
    0.5 * (idx) * ((By[i + 1][j] + By[i + 1][j - 1]) -
		   (By[i - 1][j] + By[i - 1][j - 1])) - idy * (Bx[i][j] -
							       Bx[i][j - 1]);

  // Add to energy flux

  xFlux[EN] += (by * etaPhi * jPhi - bz * etaZ * jZ);
  yFlux[EN] += (-bx * etaPhi * jPhi + bz * etaR * jR);
  zFlux[EN] += (bx * etaZ * jZ - by * etaR * jR);
}
