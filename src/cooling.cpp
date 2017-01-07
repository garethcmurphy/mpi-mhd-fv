#include "cooling.h"
#include "PhysConsts.h"

float atomic_temp_tab[50] = {
        4.04, 4.08, 4.12, 4.17, 4.21, 4.25, 4.3, 4.34, 4.38, 4.44, 4.48, 4.53,
        4.58, 4.62, 4.67, 4.71, 4.77, 4.82, 4.88, 4.93, 4.98, 5.04, 5.09, 5.14,
        5.2, 5.25, 5.3, 5.36, 5.41, 5.46, 5.51, 5.57, 5.62, 5.67, 5.73, 5.78,
        5.83, 5.88, 5.94, 5.99, 6.04, 6.09, 6.14, 6.2, 6.25, 6.3, 6.35, 6.4,
        6.45, 6.5
};


int
cooling(double *zone, double *Lcooling, double dt) {
    int hh = 0;

    double gammag = PhysConsts::gamma;
    double gammam1 = gammag - 1;
    double gammam1i = gammag - 1;
    double gammai = 1 / gammag;
    double rhoi, px, py, et, ke, p;
    double rho, velx, vely, pressure, vsnd, mach;
    double velx2;
    double vely2;
    double massfluxp;
    double massfluxn;
    double temperature;

    double ki = 24296.3696;
    double mp = 1.67262158;
    double mpi = 1.0 / mp;
    double nt = 0;

    double rate = 0, eloss, nt2, de, w, y, subdt, cv;
    double rate_temp = 0;
    double molrate = 0;
    double min_dt, lowest_temperature, subttot, elosstot;
    double e_init;
    int firststep;
    int counter;
    double nH2 = 1;
    double nh = 1;
    double chi = 1;
    double real_temp = 1;
    int rc = 0;
// Establish a system of units to work out the cooling
//     Steve's system of units
//
//     LU=10^15cm    ! YSO jet radius
//     VU=10^7 cm/s  ! 100km/s : YSO jet velocity
//     DU=10^-22     ! 120 cm^-3 100% ionised atomic H: typical ambient density
//
//     temperatureU=10^8s      ! 3.17 years
//     MU=10^23g     ! 5*10-11 solar masses
//     EU=10^37ergs  ! 26*10^3 solar luminosity integrated over 1 second
//     PU=10-8barye  ! 60 times pressuresure of reference ISM at 10^4K
//
//     proton mass= 1.67*10-47 MU
//     proton mass only used in deriving number density
//     number density only used in: 1) temperature=p/n*k 2) Lnet=n^2*LN
//     absorb 10^+47 into k and LN
//


    rho = zone[0];
    px = zone[1];
    py = zone[2];
    double pz = zone[3];
    et = zone[4];
    double bx = zone[5];
    double by = zone[6];
    double bz = zone[7];
    double b2 = bx * bx + by * by + bz * bz;

    //     cout << "cool:" ;

    // Work out temperature from zone value
    //
    rhoi = 1.0 / rho;
    velx = px * rhoi;
    vely = py * rhoi;
    double velz = pz * rhoi;
    velx2 = velx * velx;
    vely2 = vely * vely;
    double velz2 = velz * velz;
    ke = 0.5 * rho * (velx2 + vely2 + velz2);
    pressure = et - ke - 0.5 * b2;
    pressure = pressure * (gammam1);
    double pmin = 1e-5;
    if (pressure < pmin) {
        return 0;
    }
    vsnd = gammag;
    vsnd = vsnd * rhoi * pressure;
    vsnd = sqrt(vsnd);

    nt = 2 * mpi * rho;
    nt2 = nt * nt;
    temperature = ki * pressure / nt;
    real_temp = ki * pressure / nt;
    temperature = log10(temperature);
    // set up initial values

// read cooling table


// Set the minimum subcycling step
    min_dt = dt / 100.0;
    /* initialisation of cooling variables */
    lowest_temperature = atomic_temp_tab[0];    // lowest tabulated temp (LOG)
    subttot = 0.0;
    elosstot = 0.0;
    firststep = 1;

// set the initial energy
    e_init = et - ke;
    et = et - ke;

#ifdef MOLCOOL
    if (temperature < 3000)
      {
        chi = 0;
      }
    else if (temperature > 3000 && temperature < 12600)
      {
        chi = MIN ((temperature - 3000) / 7000, 0.9);
      }
    else
      {
        chi = 0;
      }
    nH2 = nt * (1 - chi);
    nh = nt * chi;
    rc = molcool (real_temp, nH2, nh, &molrate);
#endif /* MOLCOOL */
    if (temperature <= lowest_temperature) {
        rate = 0.0;
    } else {
//    need to use lookup table here
        tabfind(temperature, &rate, atomic_temp_tab);    // absorbed 10+94
    }
#ifdef MOLCOOL
    rate = rate + molrate;
#endif /* MOLCOOL */


//    initialising substep size and initial energy loss rate for all atoms
    eloss = nt2 * rate;
// loop around the substeps 
    *Lcooling = 0;


    cv = gammam1i;


    counter = 0;
    de = 0.0;
    // Substepping
    while (subttot < dt && temperature > lowest_temperature) {
        w = dt - subttot;
        y = 0.05 * et / abs(eloss);    //    max timestep for 5% change in e
        subdt = MIN (w, y);
        elosstot = eloss * subdt;
//        write(*,*) counter,temperature,'precool:t', e/elosstot
        et = et - elosstot;
        p = et / cv;

#ifdef MOLCOOL
        if (temperature < 3000)
      {
        chi = 0;
      }
        else if (temperature > 3000 && temperature < 12600)
      {
        chi = MIN ((temperature - 3000) / 7000, 0.9);
      }
        else
      {
        chi = 0;
      }
        nH2 = nt * (1 - chi);
        nh = nt * chi;
        rc = molcool (real_temp, nH2, nh, &molrate);
  //   rc = molcool (pow(10,temperature),nH2, nh  ,&rate_temp );
#endif /* MOLCOOL */
        if (temperature <= lowest_temperature) {
            eloss = 0.0;
        } else {
            tabfind(temperature, &rate, atomic_temp_tab);
#ifdef MOLCOOL
            rate = rate + molrate;
#endif /* MOLCOOL */
            eloss = nt2 * rate;
        }

        counter = counter + 1;
        subttot = subttot + subdt;

        if (firststep == 0) {
            firststep = 0;
        }

    }                //     end of substepping

//

    if (counter > 0) {
        de = et - e_init;        //     energy correction
        *Lcooling = de;
    }

    return 0;
}
