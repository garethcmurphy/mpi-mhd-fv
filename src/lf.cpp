/* $Id: riemann.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "riemann.h"
#include "PhysConsts.h"

#define DEBUG_RIEMANN
#define DEBUG_EV
#undef DEBUG_RIEMANN
//#define CONS_FLUX
#define EINFELDT_FIX
#define MINVAL(a, b) a<b?a:b
#define MAXVAL(a, b) a>b?a:b

/*
currently no entropy fix is in place 
*/
int
lf(double *leftstate,
   double *rightstate,
   double *roe_flux,
   double *res_state, int time_step, double *max_speed, int idir, int *loc) {
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
    double lstate[8];
    double rstate[8];
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

    if (rl < 0 || rr < 0) {
        cerr << "Negative density in roe solver!" << endl;
        exit(0);
    }
    outFile.open(outfilename, ofstream::app);
    if (!outFile) {
        cerr << "Unable to open file" << endl;
    }

    if (idir == 1) {
        /* Keep fluxes */
    } else if (idir == 2) {
        rl = leftstate[0];
        ul = leftstate[2];
        vl = leftstate[3];
        wl = leftstate[1];
        pl = leftstate[4];
        pl = std::max(pl, pmin);
        bul = leftstate[6];
        bvl = leftstate[7];
        bwl = leftstate[5];

        rr = rightstate[0];
        ur = rightstate[2];
        vr = rightstate[3];
        wr = rightstate[1];
        pr = rightstate[4];
        pr = std::max(pr, pmin);
        bur = rightstate[6];
        bvr = rightstate[7];
        bwr = rightstate[5];
        /* Rotate fluxes */
    } else if (idir == 3) {
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

    for (ii = 0; ii < 7; ii++) {
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
    if (fabs(cc[4]) < 1e-10)
        cc[4] = 0;
    cc[5] = (bvr - bvl);
    cc[6] = (bwr - bwl);

    /* compute the averaged quantities */
    rho_rl = 0.5 * (rr + rl);
    u_rl = 0.5 * (ul + ur);
    if (u_rl == 0) {
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
    Array2D<double> levec(7, 7);
    Array2D<double> revec(7, 7);
    Array2D<double> levc(7, 7);
    Array2D<double> revc(7, 7);

    /* Compute fast and slow speeds */
    rho = rl;
    rhoi = 1 / rho;
    srhoi = 1 / sqrt(rho);
    u = ul;
    v = vl;
    w = wl;
    bu = bul * sqrt(rhoi);
    bv = bvl * sqrt(rhoi);
    bw = bwl * sqrt(rhoi);
    bsquared = (bu * bu + bv * bv + bw * bw);
    vv2 = (u * u + v * v + w * w);
    p = pl;
    calfven2 = fabs(bu * bu);
    calfven = sqrt(calfven2);
    csound2 = gammag * p * rhoi;
    a_star2 = csound2 + bsquared;
    term = sqrt(a_star2 * a_star2 - 4.0 * csound2 * calfven2);
    term = fabs(term);
    cslow2 = 0.5 * (a_star2 - term);
    if (cfast2 < 0) {
        if (cfast2 < -1e-10) {
            cout << "Fast speed went negative!! " << cslow2 << endl;
            cout
                    << " csound2 " << csound2
                    << " calfven2 " << calfven2
                    << " bsquared " << bsquared
                    << " bx " << bu << " by " << bv << " bz " << bw << endl;
            return 1;
        } else {
            cslow2 = 0;
        }
    }
    cfast2 = 0.5 * (a_star2 + term);
    cslow = sqrt(cslow2);
    double cfast_l = sqrt(cfast2);
    if (isnan(cfast_l)) {

        std::cout << "gggg"
                  << std::endl;
        exit(0);
    }



    /* Compute fast and slow speeds */
    rho = rho_rl;
    rhoi = 1 / rho;
    srhoi = 1 / sqrt(rho);
    u = u_rl;
    v = v_rl;
    w = w_rl;
    bu = bu_rl * sqrt(rhoi);
    bv = bv_rl * sqrt(rhoi);
    bw = bw_rl * sqrt(rhoi);
    bsquared = (bu * bu + bv * bv + bw * bw);
    vv2 = (u * u + v * v + w * w);
    p = p_rl;
    calfven2 = fabs(bu * bu);
    calfven = sqrt(calfven2);
    csound2 = gammag * p * rhoi;
    a_star2 = csound2 + bsquared;
    term = sqrt(a_star2 * a_star2 - 4.0 * csound2 * calfven2);
    term = fabs(term);
    cslow2 = 0.5 * (a_star2 - term);
    cfast2 = 0.5 * (a_star2 + term);
    if (cfast2 < 0) {
        if (cfast2 < -1e-10) {
            cout << "Fast speed went negative!! " << cslow2 << endl;
            cout
                    << " csound2 " << csound2
                    << " calfven2 " << calfven2
                    << " bsquared " << bsquared
                    << " bx " << bu << " by " << bv << " bz " << bw << endl;
            return 1;
        } else {
            cslow2 = 0;
        }
    }
    cslow = sqrt(cslow2);
    cfast = sqrt(cfast2);
    if (isnan(cfast)) {
        std::cout << "gggg"
                  << std::endl;
        exit(0);
    }

    if (isinf(cfast)) {
        std::cout << "cfast inf"
                  << "p_rl  " << p_rl
                  << " csound2 " << csound2
                  << " calfven2 " << calfven2
                  << " bsquared " << bsquared
                  << " bul " << bul
                  << " bur " << bur
                  << " bu_rl " << bu_rl
                  << " bv_rl " << bv_rl
                  << " bw_rl " << bw_rl
                  << std::endl;
        exit(0);
        //cfast=0;
    }


    double clax = std::max(fabs(u_rl + cfast), fabs(u_rl - cfast));


    mass = leftstate[0];
    u = leftstate[1];
    v = leftstate[2];
    w = leftstate[3];
    p = leftstate[4];
    bu = leftstate[5];
    bv = leftstate[6];
    bw = leftstate[7];


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

    lstate[0] = mass;
    lstate[1] = mass * u;
    lstate[2] = mass * v;
    lstate[3] = mass * w;
    lstate[4] = energy;
    lstate[5] = bu;
    lstate[6] = bv;
    lstate[7] = bw;

    mass = rightstate[0];
    u = rightstate[1];
    v = rightstate[2];
    w = rightstate[3];
    p = rightstate[4];
    bu = rightstate[5];
    bv = rightstate[6];
    bw = rightstate[7];


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


    rstate[0] = mass;
    rstate[1] = mass * u;
    rstate[2] = mass * v;
    rstate[3] = mass * w;
    rstate[4] = energy;
    rstate[5] = bu;
    rstate[6] = bv;
    rstate[7] = bw;

    roe_flux[0] = 0.5 * (fl[0] + fr[0] - fabs(clax) * (rstate[0] - lstate[0]));
    roe_flux[1] = 0.5 * (fl[1] + fr[1] - fabs(clax) * (rstate[1] - lstate[1]));
    roe_flux[2] = 0.5 * (fl[2] + fr[2] - fabs(clax) * (rstate[2] - lstate[2]));
    roe_flux[3] = 0.5 * (fl[3] + fr[3] - fabs(clax) * (rstate[3] - lstate[3]));
    roe_flux[4] = 0.5 * (fl[4] + fr[4] - fabs(clax) * (rstate[4] - lstate[4]));
    roe_flux[5] = 0.;
    roe_flux[6] = 0.5 * (fl[5] + fr[5] - fabs(clax) * (rstate[6] - lstate[6]));
    roe_flux[7] = 0.5 * (fl[6] + fr[6] - fabs(clax) * (rstate[7] - lstate[7]));

    if (isnan(roe_flux[0])) {

        std::cout <<
                  " bomb "
                  << " clax= " << clax
                  << " cfast= " << cfast
                  << " u_rl= " << u_rl
                  << std::endl;
        exit(0);
    }

    return 0;
}
