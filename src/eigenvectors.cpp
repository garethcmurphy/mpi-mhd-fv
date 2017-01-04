#include "PhysConsts.h"
#include "eigenvectors.h"

/* Roe-Balsara system */

int
eigenvectors(double *sta, double **lev, double **rev, double **lec,
             double **rec, double dvy, double dvz) {
    double gammag = PhysConsts::gamma;
    int ii = 0;
    int jj = 0;
    double gammam1 = gammag - 1;
    double gammam1i = 1 / gammam1;
    double rho = 0;
    double rhoi = 0;
    double u = 0;
    double v = 0;
    double w = 0;
    double BBu = 0;
    double BBv = 0;
    double BBw = 0;
    double bu = 0;
    double bv = 0;
    double bw = 0;
    double p = 0;
    double cslow = 0;
    double calfven = 0;
    double calfven2 = 0;
    double cfast = 0;
    double cslow2 = 0;
    double cfast2 = 0;
    double csound2 = 0;
    double csound = 0;
    double term = 0;
    double a_star2 = 0;
    double bsquared = 0;
    double vv2 = 0;

    double betay = 0;
    double betaz = 0;
    double alphaf = 0;
    double alphas = 0;
    double bperp = 0;
    double icsound2 = 0;
    double icsound22 = 0;
    double phi = 0;
    double deltas = 0;
    double deltaf = 0;
    double sqrt_rho = 0;
    double sqrt_rhoi = 0;
    double pie = 3.14159;
    double vt = 0;

    double MI[7][ 7];
    double M[7][7];



    rho = sta[0];
    rhoi = 1 / rho;
    sqrt_rho = sqrt(rho);
    sqrt_rhoi = sqrt(rhoi);
    u = sta[1];
    v = sta[2];
    w = sta[3];
    BBu = sta[5];
    BBv = sta[6];
    BBw = sta[7];
    bu = sta[5] * sqrt_rhoi;
    bv = sta[6] * sqrt_rhoi;
    bw = sta[7] * sqrt_rhoi;
    p = sta[4];
    bsquared = (bu * bu + bv * bv + bw * bw);
    vv2 = (u * u + v * v + w * w);

    calfven2 = bu * bu;
    calfven = sqrt(calfven2);
    csound2 = gammag * p * rhoi;
    csound = sqrt(csound2);
    icsound2 = 1 / csound2;
    icsound22 = icsound2 * 0.5;

#ifdef DEBUG1
    cout << "eigenv: " << csound2
      << " " << gammag << " " << p << " " << rhoi << endl;
#endif /* DEBUG */

    a_star2 = (csound2 + bsquared);
    term = sqrt(a_star2 * a_star2 - 4 * csound2 * calfven2);
    cslow2 = 0.5 * (a_star2 - term);
    cfast2 = 0.5 * (a_star2 + term);
    cslow = sqrt(cslow2);
    cfast = sqrt(cfast2);

#ifdef DEBUG1
    cout << "eigenv: " << cslow2 << " " << csound2 << " " << bsquared << endl;
#endif /* DEBUG */

    /* Now must check for triple point (see Roe-Balsara 1996 */
    bperp = sqrt(bv * bv + bw * bw);
    /* Compare define a small signal speed vperp = Bperb /sqrt( rho)
     * if this is a small fraction of the sound speed then
     * */
    if ((bperp) > (0.00000001 * csound)) {
        /* Transverse field is not zero */
        betay = bv / bperp;
        betaz = bw / bperp;
    } else {
        /* Transverse field is zero   choose a random direction - Brio and Wu */
        /* These quantities indeterminate if transverse field = 0
         *
         * */
        vt = sqrt(dvy * dvy + dvz * dvz);
        if (vt > (0.00000001 * csound)) {
            betay = sgn(bu) * dvy / vt;
            betaz = sgn(bu) * dvz / vt;
        } else {
            betay = 1 / sqrt(2.0);
            betaz = 1 / sqrt(2.0);
            betay = 0;
            betaz = 0;
        }
    }

    if (cfast > (cslow + 0.0000000001 * csound)) {
//      Not degenerate case
        alphaf = (csound2 - cslow2) / (cfast2 - cslow2);
        alphas = (cfast2 - csound2) / (cfast2 - cslow2);
        alphaf = sqrt(abs(alphaf));
        alphas = sqrt(abs(alphas));
    } else if (fabs(calfven - csound) > 1e-7 * csound) {
        phi = atan(bperp / (fabs(bu) - csound));
        alphas = fabs(cos(phi / 2)) + deltas;
        alphaf = fabs(sin(phi / 2)) + deltaf;
    } else {
        phi = 0.25 * pie * sgn(calfven - csound);
        alphas = fabs(cos(phi));
        alphaf = fabs(sin(phi));

    }

/*
	 if ( cfast -cslow == 0) {
		 alphaf=1.0;
		 alphas=0.0;
	 } else if ( csound - cslow <=0.){
		 alphaf=1.0;
		 alphas=0.0;
	 } else if ( cslow - csound <=0.){
		 alphaf=1.0;
		 alphas=0.0;
	 } else {
      alphaf = (csound2 - cslow2) / (cfast2 - cslow2);
      alphas = (cfast2 - csound2) / (cfast2 - cslow2);
      alphaf = sqrt (abs (alphaf));
      alphas = sqrt (abs (alphas));
	 }
*/


#ifdef DEBUG1
    cout << "eigenv: " << dels
      << " " << delf << " " << betas << " " << betaf << endl;
    exit (0);
#endif /* DEBUG */

#ifdef DEBUG1
    cout << "eigenv: " << dels
      << " " << rho << " " << cslow2 << " " << BBu << endl;
#endif /* DEBUG */

    /* compute right eigenvectors */
    /* r u-cfast */
    rev[0][0] = alphaf * rho;
    rev[1][0] = (-alphaf * cfast);
    rev[2][0] = (alphas * cslow * betay * sgn(bu));
    rev[3][0] = (alphas * cslow * betaz * sgn(bu));
    rev[4][0] = alphaf * rho * csound2;
    rev[5][0] = alphas * sqrt_rho * csound * betay;
    rev[6][0] = alphas * sqrt_rho * csound * betaz;
    /* r u-calfven */
    rev[0][1] = 0;
    rev[1][1] = 0;
    rev[2][1] = (-betaz);
    rev[3][1] = (betay);
    rev[4][1] = 0;
    rev[5][1] = (-sqrt_rho) * betaz * sgn(bu);
    rev[6][1] = sqrt_rho * betay * sgn(bu);
    /* r u-cslow */
    rev[0][2] = alphas * rho;
    rev[1][2] = (-alphas * cslow);
    rev[2][2] = (-alphaf * cfast * betay * sgn(bu));
    rev[3][2] = (-alphaf * cfast * betaz * sgn(bu));
    rev[4][2] = alphas * rho * csound2;
    rev[5][2] = (-alphaf) * sqrt_rho * csound * betay;
    rev[6][2] = (-alphaf) * sqrt_rho * csound * betaz;
    /* r u */
    rev[0][3] = 1;
    rev[1][3] = 0;
    rev[2][3] = 0;
    rev[3][3] = 0;
    rev[4][3] = 0;
    rev[5][3] = 0;
    rev[6][3] = 0;
    /* r u+cslow */
    rev[0][4] = alphas * rho;
    rev[1][4] = alphas * cslow;
    rev[2][4] = (alphaf * cfast * betay * sgn(bu));
    rev[3][4] = (alphaf * cfast * betaz * sgn(bu));
    rev[4][4] = alphas * rho * csound2;
    rev[5][4] = (-alphaf) * sqrt_rho * csound * betay;
    rev[6][4] = (-alphaf) * sqrt_rho * csound * betaz;
    /* r u+calfven */
    rev[0][5] = 0;
    rev[1][5] = 0;
    rev[2][5] = betaz;
    rev[3][5] = (-betay);
    rev[4][5] = 0;
    rev[5][5] = (-sqrt_rho) * betaz * sgn(bu);
    rev[6][5] = sqrt_rho * betay * sgn(bu);
    /* r u+cfast */
    rev[0][6] = alphaf * rho;
    rev[1][6] = alphaf * cfast;
    rev[2][6] = -(alphas * cslow * betay * sgn(bu));
    rev[3][6] = -(alphas * cslow * betaz * sgn(bu));
    rev[4][6] = alphaf * rho * csound2;
    rev[5][6] = alphas * sqrt_rho * csound * betay;
    rev[6][6] = alphas * sqrt_rho * csound * betaz;


    /* compute left eigenvectors */
    /* r u-cfast */
    lev[0][0] = 0;
    lev[1][0] = icsound22 * (-alphaf * cfast);
    lev[2][0] = icsound22 * (alphas * cslow * betay * sgn(bu));
    lev[3][0] = icsound22 * (alphas * cslow * betaz * sgn(bu));
    lev[4][0] = icsound22 * alphaf * rhoi;
    lev[5][0] = icsound22 * alphas * sqrt_rhoi * csound * betay;
    lev[6][0] = icsound22 * alphas * sqrt_rhoi * csound * betaz;
    /* r u-calfven */
    lev[0][1] = 0;
    lev[1][1] = 0;
    lev[2][1] = 0.5 * (-betaz);
    lev[3][1] = 0.5 * (betay);
    lev[4][1] = 0;
    lev[5][1] = 0.5 * (-betaz * sgn(bu) * sqrt_rhoi);
    lev[6][1] = 0.5 * (betay * sgn(bu) * sqrt_rhoi);
    /* r u-cslow */
    lev[0][2] = 0;
    lev[1][2] = icsound22 * (-alphas * cslow);
    lev[2][2] = icsound22 * (-alphaf * cfast * betay * sgn(bu));
    lev[3][2] = icsound22 * (-alphaf * cfast * betaz * sgn(bu));
    lev[4][2] = icsound22 * alphas * rhoi;
    lev[5][2] = icsound22 * (-alphaf) * sqrt_rhoi * csound * betay;
    lev[6][2] = icsound22 * (-alphaf) * sqrt_rhoi * csound * betaz;
    /* r u */
    lev[0][3] = 1;
    lev[1][3] = 0;
    lev[2][3] = 0;
    lev[3][3] = 0;
    lev[4][3] = -1 / csound2;
    lev[5][3] = 0;
    lev[6][3] = 0;
    /* r u+cslow */
    lev[0][4] = 0;
    lev[1][4] = icsound22 * alphas * cslow;
    lev[2][4] = icsound22 * (alphaf * cfast * betay * sgn(bu));
    lev[3][4] = icsound22 * (alphaf * cfast * betaz * sgn(bu));
    lev[4][4] = icsound22 * alphas * rhoi;
    lev[5][4] = icsound22 * (-alphaf) * sqrt_rhoi * csound * betay;
    lev[6][4] = icsound22 * (-alphaf) * sqrt_rhoi * csound * betaz;
    /* r u+calfven */
    lev[0][5] = 0;
    lev[1][5] = 0;
    lev[2][5] = 0.5 * betaz;
    lev[3][5] = 0.5 * (-betay);
    lev[4][5] = 0;
    lev[5][5] = 0.5 * (-betaz * sgn(bu) * sqrt_rhoi);
    lev[6][5] = 0.5 * (betay * sgn(bu) * sqrt_rhoi);
    /* r u+cfast */
    lev[0][6] = 0;
    lev[1][6] = icsound22 * alphaf * cfast;
    lev[2][6] = icsound22 * (-alphas * cslow * betay * sgn(bu));
    lev[3][6] = icsound22 * (-alphas * cslow * betaz * sgn(bu));
    lev[4][6] = icsound22 * alphaf * rhoi;
    lev[5][6] = icsound22 * alphas * sqrt_rhoi * csound * betay;
    lev[6][6] = icsound22 * alphas * sqrt_rhoi * csound * betaz;

    for (ii = 0; ii < 7; ii++) {
        for (jj = 0; jj < 7; jj++) {
            if (isnan(lev[ii][jj])) {
                cout << "lev[" << ii << "," << jj << "] is nan " << endl;
                exit(0);
				                if (isnan(rev[ii][jj])) {
                    cout << "rev[" << ii << "," << jj << "] is nan " << endl;
                    exit(0);
            }
        }
    }
}


//-----Transform primitive eigenvectors to conservative eigen vectors;
//     See Hirsch Vol. II p.147 eqns.[16.2.47];
//;

M[0][0] = 1.0;
//       M[1][2] = 0.0;
//       M[1][3] = 0.0;
//       M[1][4] = 0.0;
//       M[1][5] = 0.0;
//       M[1][6] = 0.0;
//       M[1][7] = 0.0;


M[1][0] =
u;
M[1][1] =
rho;
//       M[2][3] = 0.0;
//       M[2][4] = 0.0;
//       M[2][5] = 0.0;
//       M[2][6] = 0.0;
//       M[2][7] = 0.0;


M[2][0] =
v;
//       M[3][2] = 0.0;
M[2][2] =
rho;
//       M[3][4] = 0.0;
//       M[3][5] = 0.0;
//       M[3][6] = 0.0;
//       M[3][7] = 0.0;


M[3][0] =
w;
//       M[4][2] = 0.0;
//       M[4][3] = 0.0;
M[3][3] =
rho;
//       M[4][5] = 0.0;
//       M[4][6] = 0.0;
//       M[4][7] = 0.0;


M[4][0] = 0.5 *
vv2;
M[4][1] =
rho *u;
M[4][2] =
rho *v;
M[4][3] =
rho *w;
M[4][4] =
gammam1i;
M[4][5] =
BBv;
M[4][6] =
BBw;


//       M[6][1] = 0.0;
//       M[6][2] = 0.0;
//       M[6][3] = 0.0;
//       M[6][4] = 0.0;
//       M[6][5] = 0.0;
M[5][5] = 1.0;
//       M[6][7] = 0.0;


//       M[7][1] = 0.0;
//       M[7][2] = 0.0;
//       M[7][3] = 0.0;
//       M[7][4] = 0.0;
//       M[7][5] = 0.0;
//       M[7][6] = 0.0;
M[6][6] = 1.0;
//;

MI[0][0] = 1.0;
MI[1][0] = -
rhoi *u;
MI[2][0] = -
rhoi *v;
MI[3][0] = -
rhoi *w;

MI[4][0] = 0.5 *
gammam1 *vv2;
//       MI[6][1] = 0.0;
//       MI[7][1] = 0.0;


//       MI[1][2] = 0.0;
MI[1][1] =
rhoi;
//       MI[3][2] = 0.0;
//       MI[4][2] = 0.0;
MI[4][1] = -
gammam1 *u;
//       MI[6][2] = 0.0;
//       MI[7][2] = 0.0;


//       MI[1][3] = 0.0;
//       MI[2][3] = 0.0;
MI[2][2] =
rhoi;
//       MI[4][3] = 0.0;
MI[4][2] = -
gammam1 *v;
//       MI[6][3] = 0.0;
//       MI[7][3] = 0.0;


//       MI[1][4] = 0.0;
//       MI[2][4] = 0.0;
//       MI[3][4] = 0.0;
MI[3][3] =
rhoi;
MI[4][3] = -
gammam1 *w;
//       MI[6][4] = 0.0;
//       MI[7][4] = 0.0;


//       MI[1][5] = 0.0;
//       MI[2][5] = 0.0;
//       MI[3][5] = 0.0;
//       MI[4][5] = 0.0;
MI[4][4] =
gammam1;
//       MI[6][5] = 0.0;
//       MI[7][5] = 0.0;


//       MI[1][6] = 0.0;
//       MI[2][6] = 0.0;
//       MI[3][6] = 0.0;
//       MI[4][6] = 0.0;
MI[4][5] = -
gammam1 *BBv;
MI[5][5] = 1.0;
//       MI[7][6] = 0.0;


//       MI[1][7] = 0.0;
//       MI[2][7] = 0.0;
//       MI[3][7] = 0.0;
//       MI[4][7] = 0.0;
MI[4][6] = -
gammam1 *BBw;
//       MI[6][7] = 0.0;
MI[6][6] = 1.0;
//;
//-----compute conserved eigenvectors;
//;
//
//
int k=0;
for (k = 0; k < 7; k++) 
{
rec[0][k] = M[0][0] * rev[0][k];
rec[1][k] = M[1][0] * rev[0][k] + M[1][1] * rev[1][k];
rec[2][k] = M[2][1] * rev[0][k] + M[2][2] * rev[2][k];
rec[3][k] = M[3][0] * rev[0][k] + M[3][3] * rev[3][k];
rec[4][k] = M[4][0] * rev[0][k] + M[4][1] * rev[1][k]
+ M[4][2] * rev[2][k] + M[4][3] * rev[3][k]
+ M[4][4] * rev[4][k] + M[4][5] * rev[5][k] + M[4][6] * rev[6][k];
rec[5][k] = M[5][5] * rev[5][k];
rec[6][k] = M[6][6] * rev[6][k];


/*
   lec[k][0] = lev[k][0] * MI[0][0] + lev[k][1] * MI[1][0] + lev[k][2] * MI[2][0] + lev[k][3] * MI[3][0] + lev[k][4] * MI[4][0];
   lec[k][1] = lev[k][1] * MI[1][1] + lev[k][4] * MI[4][1];
   lec[k][2] = lev[k][2] * MI[2][2] + lev[k][4] * MI[4][2];
   lec[k][3] = lev[k][3] * MI[3][3] + lev[k][4] * MI[4][3];
   lec[k][4] = lev[k][4] * MI[4][4];
   lec[k][5] = lev[k][4] * MI[4][5] + lev[k][5] * MI[5][5];
   lec[k][6] = lev[k][4] * MI[4][6] + lev[k][6] * MI[6][6];

   lec[0][k] = lev[k][0] * MI[0][0] + lev[k][1] * MI[1][0] + lev[k][2] * MI[2][0] + lev[k][3] * MI[3][0] + lev[k][4] * MI[4][0];
   lec[1][k] = lev[k][1] * MI[1][1] + lev[k][4] * MI[4][1];
   lec[2][k] = lev[k][2] * MI[2][2] + lev[k][4] * MI[4][2];
   lec[3][k] = lev[k][3] * MI[3][3] + lev[k][4] * MI[4][3];
   lec[4][k] = lev[k][4] * MI[4][4];
   lec[5][k] = lev[k][4] * MI[4][5] + lev[k][5] * MI[5][5];
   lec[6][k] = lev[k][4] * MI[4][6] + lev[k][6] * MI[6][6];
 */


lec[0][k] =
lev[0][k] * MI[0][0] + lev[1][k] * MI[1][0] + lev[2][k] * MI[2][0] +
lev[3][k] * MI[3][0] + lev[4][k] * MI[4][0];
lec[1][k] = lev[1][k] * MI[1][1] + lev[4][k] * MI[4][1];
lec[2][k] = lev[2][k] * MI[2][2] + lev[4][k] * MI[4][2];
lec[3][k] = lev[3][k] * MI[3][3] + lev[4][k] * MI[4][3];
lec[4][k] = lev[4][k] * MI[4][4];
lec[5][k] = lev[4][k] * MI[4][5] + lev[5][k] * MI[5][5];
lec[6][k] = lev[4][k] * MI[4][6] + lev[6][k] * MI[6][6];
}


return 0;
}

int
eigenvalues(double *sta, double *eigenval) {
    double gammag = PhysConsts::gamma;
    int ii = 0;
    int jj = 0;
    double gammam1 = gammag - 1;
    double gammam1i = 1 / gammam1;
    double rho = 0;
    double rhoi = 0;
    double u = 0;
    double v = 0;
    double w = 0;
    double BBu = 0;
    double BBv = 0;
    double BBw = 0;
    double bu = 0;
    double bv = 0;
    double bw = 0;
    double p = 0;
    double cslow = 0;
    double calfven = 0;
    double calfven2 = 0;
    double cfast = 0;
    double cslow2 = 0;
    double cfast2 = 0;
    double csound2 = 0;
    double csound = 0;
    double term = 0;
    double a_star2 = 0;
    double bsquared = 0;
    double vv2 = 0;
    double dels = 0;
    double delf = 0;
    double betaf = 0;
    double betas = 0;
    double sqrt2 = sqrt(2.0);

    double betay = 0;
    double betaz = 0;
    double alphaf = 0;
    double alphas = 0;
    double bperp = 0;
    double sqrt_rho = 0;
    double sqrt_rhoi = 0;
	double icsound2 =0;
	double icsound22 =0;
		

    int k = 0;
    rho = sta[0];
    rhoi = 1 / rho;
    sqrt_rho = sqrt(rho);
    sqrt_rhoi = sqrt(rhoi);
    u = sta[1];
    v = sta[2];
    w = sta[3];
    BBu = sta[5];
    BBv = sta[6];
    BBw = sta[7];
    bu = sta[5] * sqrt_rhoi;
    bv = sta[6] * sqrt_rhoi;
    bw = sta[7] * sqrt_rhoi;
    p = sta[4];
    bsquared = (bu * bu + bv * bv + bw * bw);
    vv2 = (u * u + v * v + w * w);

    calfven2 = bu * bu;
    calfven = sqrt(calfven2);
    csound2 = gammag * p * rhoi;
    csound = sqrt(csound2);
    icsound2 = 1 / csound2;
    icsound22 = icsound2 * 0.5;

#ifdef DEBUG1
    cout << "eigenv: " << csound2
      << " " << gammag << " " << p << " " << rhoi << endl;
#endif /* DEBUG */

    a_star2 = (csound2 + bsquared);
    term = sqrt(a_star2 * a_star2 - 4 * csound2 * calfven2);
    cslow2 = 0.5 * (a_star2 - term);
    cfast2 = 0.5 * (a_star2 + term);
    cslow = sqrt(cslow2);
    cfast = sqrt(cfast2);
    eigenval[0] = u - cfast;
    eigenval[1] = u - calfven;
    eigenval[2] = u - cslow;
    eigenval[3] = u;
    eigenval[4] = u + cslow;
    eigenval[5] = u + calfven;
    eigenval[6] = u + cfast;
    return 0;
}
