
#include "parupdate.h"
#include "physics.h"
#include "cooling.h"

double min4(
        double a,
        double b,
        double c,
        double d
) {
    return std::min(std::min(a, b), std::min(c, d));

}

double
getpres(double *array) {


    double rho, rhoi;
    double px;
    double py, pz, et, velx, vely, velz, velx2, vely2, velz2, ke;
    double bx, by, bz, bsquared;
    double gammag = PhysConsts::gamma;
    double gammam1 = gammag - 1;
    double gammam1i = 1.0 / gammam1;

    rho = array[0];
    px = array[1];
    py = array[2];
    pz = array[3];
    et = array[4];
    bx = array[5];
    by = array[6];
    bz = array[7];
    rhoi = 1.0 / rho;
    velx = px * rhoi;
    vely = py * rhoi;
    velz = pz * rhoi;
    velx2 = velx * velx;
    vely2 = vely * vely;
    velz2 = velz * velz;
    ke = 0.5 * rho * (velx2 + vely2 + velz2);
    bsquared = bx * bx + by * by + bz * bz;
    double pressure = et - ke - 0.5 * bsquared;
    pressure = pressure * (gammam1);

    return 0;
}


int
parupdate(TNT::Array2D<unk> newmesh,
          TNT::Array2D<unk> sourcemesh,
          TNT::Array2D<unk> mesh,
          TNT::Array2D<flux> xflux,
          TNT::Array2D<flux> yflux,
          TNT::Array2D<flux> SecondOrdCorr,
          TNT::Array2D<double> faceBx,
          TNT::Array2D<double> faceBy,
          TNT::Array2D<double> xjx,
          TNT::Array2D<double> xjy,
          TNT::Array2D<double> xjz,
          TNT::Array2D<double> yjx,
          TNT::Array2D<double> yjy,
          TNT::Array2D<double> yjz,
          TNT::Array2D<double> faceEx,
          TNT::Array2D<double> faceEy,
          TNT::Array1D<double> SoundSpeedMidplane,
          double delta_t, double delta_x, double delta_y, MPI_Comm Cart_comm,
          int second) {
    TNT::Array2D<double> oldfaceBx;
    TNT::Array2D<double> oldfaceBy;
    oldfaceBx = faceBx.copy();
    oldfaceBy = faceBy.copy();
//      double dtody=delta_x* dtodx /delta_y ;
//      double delta_t =delta_x* dtodx  ; 
    double dtodx = delta_t / delta_x;
    double dtody = delta_t / delta_y;
    double idx = 1.0 / delta_x;
    double idy = 1.0 / delta_y;

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    double pmin = 1e-18;
    int nx = mesh.dim1();
    int ny = mesh.dim2();
    TNT::Array2D<double> pressurearr(nx, ny);


    int myaddress[] = {0, 0};

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    int ndims = 2;
    int coords[2];
    MPI_Cartdim_get(Cart_comm, &ndims);
    MPI_Cart_coords(Cart_comm, myid, ndims, coords);
    myaddress[0] = (nx - 2 * NGC) * coords[0];
    myaddress[1] = (ny - 2 * NGC) * coords[1];

    // Sweep for time update
    for (int ii = 2; ii < nx - 2; ii++) {
        for (int jj = 2; jj < ny - 2; jj++) {

            if (isnan(mesh[ii][jj]_MASS)) {
                std::cout
                        << " proc: " << myid
                        << " Mass is nan"
                        << "(" << ii << "," << jj << ")" << std::endl;
                std::cout
                        << __FUNCTION__
                        << ": "
                        << " " << mesh[ii][jj + 1]_MASS
                        << " " << mesh[ii][jj - 1]_MASS
                        << " " << mesh[ii + 1][jj]_MASS
                        << " " << mesh[ii - 1][jj]_MASS
                        << " " << dtodx
                        << " " << yflux[ii][jj + 1]_MASS
                        << " " << yflux[ii][jj - 1]_MASS
                        << " " << xflux[ii + 1][jj]_MASS
                        << " " << xflux[ii - 1][jj]_MASS << std::endl;
                MPI_Finalize();
                exit(0);
            }


            for (int hh = 0; hh < NE; hh++) {
#ifdef CYLINDRICAL
                double myx = (double) myaddress[0] + (ii - NGC);

                double lm = myx;
                double lp = lm + 1;
                double dri = 0;

                lm = (double) myaddress[0] + ii - NGC;
                lp = lm + 1.0;
                dri = 2.0 / (2.0 * lp + 1);

                newmesh[ii][jj].array[hh] =
              mesh[ii][jj].array[hh]
              - dtodx * ((lp) * xflux[ii + 1][jj].array[hh] -
                     (lm) * xflux[ii][jj].array[hh]) * dri -
              dtody * (yflux[ii][jj + 1].array[hh] -
                   yflux[ii][jj].array[hh]);
#else


                /*
                if (second)
                {

              newmesh[ii][jj].array[hh] =
            0.5*(mesh[ii][jj].array[hh]
                    +sourcemesh[ii][jj].array[hh]
                    )
            - dtodx * (xflux[ii + 1][jj].array[hh] -
                   xflux[ii][jj].array[hh]) - dtody * (yflux[ii][jj +
                                         1].
                                       array[hh] -
                                       yflux[ii][jj].
                                       array[hh]);
                }
                else
                    */
                {

                    newmesh[ii][jj].array[hh] =
                            mesh[ii][jj].array[hh]
                            - dtodx * (xflux[ii + 1][jj].array[hh] - xflux[ii][jj].array[hh])
                            - dtody * (yflux[ii][jj + 1].array[hh] - yflux[ii][jj].array[hh]);
                }

#endif
            }


//#define CHECKSYM
#ifdef CHECKSYM
            if (
                (ii==11 && jj==11)
                )
            {
            int 	hh=0;

            std::cout  << std::setprecision(10)
            << " "  <<
          newmesh[ii][jj].array[hh]
            << " "  <<
        mesh[ii][jj].array[hh]
            << " "  <<
        - dtodx
            << " "  <<
        xflux[ii + 1][jj].array[hh]
            << " "  <<
        xflux[ii][jj].array[hh]
            << " "  <<
         dtody
            << " "  <<
         yflux[ii][jj + 1].  array[hh]
            << " "  <<
         yflux[ii][jj].  array[hh]
         << std::endl;
            }
#endif

#ifdef CHECK_NEG_DEN
            bool NegDensity = newmesh[ii][jj] _MASS < 0 ? true : false;
            bool NegEnergy = newmesh[ii][jj] _ENER < 0 ? true : false;

            if (NegDensity || NegEnergy)
              {
                int myid;
                MPI_Comm_rank (MPI_COMM_WORLD, &myid);

                if (NegDensity)
              {
                std::cout
                  << " proc: " << myid
                  << " density is < 0 "
                  << "(" << ii << "," << jj << ")"
                  << " density =  " << newmesh[ii][jj] _MASS << std::endl;
              }
                else
              {
                std::cout
                  << " proc: " << myid
                  << " energy is < 0 "
                  << "(" << ii << "," << jj << ")"
                  << " energy =  " << newmesh[ii][jj] _ENER << std::endl;
              }


                std::cout << __FUNCTION__ << ": " << std::endl;
                std::cout << "i+1,j " << std::endl;
                for (int hh = 0; hh < 8; ++hh)
              {
                std::cout << " " << mesh[ii + 1][jj].array[hh];
              }
                std::cout << std::endl;
      //
                std::cout << "i-1,j " << std::endl;
                for (int hh = 0; hh < 8; ++hh)
              {
                std::cout << " " << mesh[ii - 1][jj].array[hh];
              }
                std::cout << std::endl;
      //
                std::cout << "i,j+1 " << std::endl;
                for (int hh = 0; hh < 8; ++hh)
              {
                std::cout << " " << mesh[ii][jj + 1].array[hh];
              }
                std::cout << std::endl;
      //
                std::cout << "i,j-1 " << std::endl;
                for (int hh = 0; hh < 8; ++hh)
              {
                std::cout << " " << mesh[ii][jj - 1].array[hh];
              }
                std::cout << std::endl;
      //
                std::cout << "update " << std::endl;
                for (int hh = 0; hh < 8; ++hh)
              {
                std::cout
                  << " " << newmesh[ii][jj].array[hh]
                  << " " << mesh[ii][jj].array[hh]
                  << " " << dtodx
                  << " " << xflux[ii + 1][jj].array[hh]
                  << " " << xflux[ii][jj].array[hh]
                  << " " << dtody
                  << " " << yflux[ii][jj + 1].array[hh]
                  << " " << yflux[ii][jj].array[hh] << std::endl;
              }
                MPI_Abort (MPI_COMM_WORLD, 5678);
                exit (0);
              }
#endif


            double rho, rhoi;
            double px;
            double py, pz, et, velx, vely, velz, velx2, vely2, velz2, ke;    //, pressure;
            double bx, by, bz, bsquared;
            double gammag = PhysConsts::gamma;
            double gammam1 = gammag - 1;
            double gammam1i = 1.0 / gammam1;

            rho = newmesh[ii][jj].array[0];
            px = newmesh[ii][jj].array[1];
            py = newmesh[ii][jj].array[2];
            pz = newmesh[ii][jj].array[3];
            et = newmesh[ii][jj].array[4];
            bx = newmesh[ii][jj].array[5];
            by = newmesh[ii][jj].array[6];
            bz = newmesh[ii][jj].array[7];
            rhoi = 1.0 / rho;
            velx = px * rhoi;
            vely = py * rhoi;
            velz = pz * rhoi;
            velx2 = velx * velx;
            vely2 = vely * vely;
            velz2 = velz * velz;
            ke = 0.5 * rho * (velx2 + vely2 + velz2);
            bsquared = bx * bx + by * by + bz * bz;
            pressurearr[ii][jj] = et - ke - 0.5 * bsquared;
            pressurearr[ii][jj] = pressurearr[ii][jj] * (gammam1);
#ifdef INJECT_THERMAL_ENERGY
            if (pressurearr[ii][jj] < 0)
              {



                newmesh[ii][jj].array[4] = pmin * gammam1i + bsquared + ke;
              }
#endif
            if (0 && pressurearr[ii][jj] < 0) {

                int myid;
                MPI_Comm_rank(MPI_COMM_WORLD, &myid);
                std::cout
                        << " proc: " << myid
                        << " pressure is < 0 "
                        << "(" << ii << "," << jj << ")"
                        << " pressure =  " << pressurearr[ii][jj] << std::endl;
                std::cout << __FUNCTION__ << ": " << std::endl;
                for (int hh = 0; hh < 8; ++hh) {
                    std::cout << " " << mesh[ii + 1][jj].array[hh];
                }
                std::cout << std::endl;
//
                for (int hh = 0; hh < 8; ++hh) {
                    std::cout << " " << mesh[ii - 1][jj].array[hh];
                }
                std::cout << std::endl;
//
                for (int hh = 0; hh < 8; ++hh) {
                    std::cout << " " << mesh[ii][jj + 1].array[hh];
                }
                std::cout << std::endl;
//
                for (int hh = 0; hh < 8; ++hh) {
                    std::cout << " " << mesh[ii][jj - 1].array[hh];
                }
                std::cout << std::endl;
//
                for (int hh = 0; hh < 8; ++hh) {
                    std::cout
                            << " " << newmesh[ii][jj].array[hh]
                            << " " << mesh[ii][jj].array[hh]
                            << " " << dtodx
                            << " " << yflux[ii][jj + 1].array[hh]
                            << " " << yflux[ii][jj].array[hh]
                            << " " << xflux[ii + 1][jj].array[hh]
                            << " " << xflux[ii][jj].array[hh] << std::endl;
                }
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Abort(MPI_COMM_WORLD, 5678);
                exit(0);
//          
            }

        }
    }


// Source terms to account for non-conservation of radial momentum
#ifdef CYLINDRICAL


    for (int ii = 2; ii < nx - 2; ii++)
      {
        for (int jj = 2; jj < ny - 2; jj++)
      {
        // What is location
        double loc[2];
        loc[0] = (ii - NGC) + myaddress[0];
        loc[1] = (jj - NGC) + myaddress[1];
        double gammam1 = PhysConsts::gamma - 1;
        double rr = sourcemesh[ii][jj] _MASS;
        double px = sourcemesh[ii][jj] _MOMX;
        double py = sourcemesh[ii][jj] _MOMY;
        double pz = sourcemesh[ii][jj] _MOMZ;
        double et = sourcemesh[ii][jj] _ENER;
        double bx = sourcemesh[ii][jj] _B_X;
        double by = sourcemesh[ii][jj] _B_Y;
        double bz = sourcemesh[ii][jj] _B_Z;
        double ri = 1.0 / rr;
        double vx = px * ri;
        double vy = py * ri;
        double vz = pz * ri;
        double ke = 0.5 * rr * (vx * vx + vy * vy + vz * vz);
        double b2 = 0.5 * (bx * bx + by * by + bz * bz);
        double pressure = et - ke - b2;
        double ptot = et - ke;
        pressure = pressure * gammam1;
        pressure = std::max (pmin, pressure);
        double lm = loc[0];
        double lp = loc[0] + 1;
        double dri = 1 / ((loc[0] + 0.5));
        lm = (double) myaddress[0] + ii - 2.0;
        lp = lm + 1.0;
        dri = 2.0 / (2.0 * lp + 1.0);
        double rgi = (lm*lm + lm + 1.0 /3.0) * delta_x / (lm +0.5);
        //rgi = (lm * lm - lm + 1.0 / 3.0) * delta_x / (lm - 0.5);

        if (ptot < 0)
          {
            std::cout
              << " " << ii <<  " " << jj
              << " ptot <0" << std::endl;
          }
        //if (second == 0)
        if (1 )
          {
  //#define NOSOURCE
#ifndef NOSOURCE
  /*
   * Geometrical source terms from
   * NRL formulary page 8 + MHD Equations
   *
   * */
            newmesh[ii][jj] _MOMX +=
          dtodx * (rr * vz * vz + pressure + b2 - bz * bz) * dri;
          //dtodx * (rr * vz * vz + 0 - bz * bz) * dri;

            newmesh[ii][jj] _MOMZ +=
          dtodx * (-rr * vx * vz + 0 + bx * bz) * dri;

            newmesh[ii][jj] _B_Z += dtodx * (vx * bz - vz * bx) * dri;
#endif
          }
        else
          {
            double dmdr = SecondOrdCorr[ii][jj] _MASS;
            double dvxdr = SecondOrdCorr[ii][jj] _MOMX;
            double dvydr = SecondOrdCorr[ii][jj] _MOMY;
            double dvzdr = SecondOrdCorr[ii][jj] _MOMZ;
            double dpdr = SecondOrdCorr[ii][jj] _ENER;
            double dbxdr = SecondOrdCorr[ii][jj] _B_X;
            double dbydr = SecondOrdCorr[ii][jj] _B_Y;
            double dbzdr = SecondOrdCorr[ii][jj] _B_Z;

            double svx =
          dmdr * vz * vz + 2 * rr * dvzdr * vz + dpdr + bx * dbxdr +
          by * dbydr - bz * dbzdr;
            double svz =
          dmdr * vx * vz + rr * dvxdr * vz + rr * vx * dvzdr -
          bx * dbzdr - bz * dbxdr;
            double sbz = vx * dbzdr + bz * dvxdr - vz * dbxdr - bx * dvzdr;
            newmesh[ii][jj] _MOMX +=
          dtodx * (rr * vz * vz + pressure + b2 - bz * bz +
               (delta_x * (ii - 0.5) - rgi) * svx) * dri;
            newmesh[ii][jj] _MOMZ -=
          dtodx * (rr * vx * vz - bx * bz +
               (delta_x * (ii - 0.5) - rgi) * svz) * dri;
            newmesh[ii][jj] _B_Z +=
          dtodx * (vx * bz - vz * bx +
               (delta_x * (ii - 0.5) - rgi) * sbz) * dri;
          }
      }
      }

#else
    if (myid == 0) { std::cout << "No geometrical source term." << std::endl; }
#endif





/* 
 *
 * Flux constrained transport
 *
 */

    TNT::Array2D<double> Ez(nx, ny);
    Ez = 9.e99;
    //Ez = 0;
    // Determine Ez at cell corners
    // Average face-centred electric field onto cell corners
    // Use upwinding?
    // Some of the values of Ez are based on dodgy values of faceEx and faceEy
    // i.e. the boundaries are not good values
    // fluxes are determined from 1,1 to nx-2, ny -2
    // same for faceEx and faceEy
    // but here we use face Ex and faceEy from 0,1 and 1,0
    // copy faceey and faceEx values in using emf_exchange
#ifdef MAES
#ifndef  RESISTIVE
    if (myid == 0) { std::cout << "No resisitivity in MAES" << std::endl; }
#endif
#endif
    /*
#ifdef PERIODIC
 for (int i = 0; i < nx ; ++i)
 {
     faceEx[i][ny-2]=faceEx[i][2] ;
     faceEx[i][1]= (faceEx[i][ny-3]) ;

 }
 for (int j = 1; j < ny - 1; ++j)
 {
     faceEy[nx-2][j]= faceEy[2][j];
     faceEy[1][j]= (faceEy[nx-3][j]) ;
 }
#else
 for (int i = 0; i < nx ; ++i)
 {
     faceEx[i][ny-2]=faceEx[i][ny-3] ;
     faceEx[i][1]= (faceEx[i][2]) ;
 }
 for (int j = 1; j < ny - 1; ++j)
 {
     faceEy[nx-2][j]= faceEy[nx-3][j];
     faceEy[1][j]= (-faceEy[2][j]) ;
 }
#endif
 */



    for (int i = 2; i < nx - 1; ++i)
        for (int j = 2; j < ny - 1; ++j) {

            // Balsara and Spicer style upwinding
            // See Steve O'Sullivan's thesis
            double alphaB = 0.5;
            bool swa = false;
            bool swb = 0;

            double dp1 = pressurearr[i + 1][j + 1] - pressurearr[i][j + 1];
            double dp2 = pressurearr[i + 1][j] - pressurearr[i][j];
            double dp3 = pressurearr[i + 1][j + 1] - pressurearr[i + 1][j];
            double dp4 = pressurearr[i + 1][j] - pressurearr[i + 1][j - 1];
            double delpz = 0.5 * fabs(dp1 - dp2);
            double delpr = 0.5 * fabs(dp3 - dp4);

            if (delpz + delpr > alphaB * min4(pressurearr[i][j], pressurearr[i + 1][j], pressurearr[i][j + 1],
                                              pressurearr[i + 1][j + 1])) {
                swa = true;
            }
            double betaB = 0.5;

            double rho1i = 1.0 / mesh[i][j]_MASS;
            double rho2i = 1.0 / mesh[i][j]_MASS;
            double rho3i = 1.0 / mesh[i][j]_MASS;
            double rho4i = 1.0 / mesh[i][j]_MASS;

            double b21 =
                    mesh[i][j]_B_X * mesh[i][j]_B_X + mesh[i][j]_B_Y * mesh[i][j]_B_Y + mesh[i][j]_B_Z * mesh[i][j]_B_Z;
            double b22 = mesh[i + 1][j]_B_X * mesh[i + 1][j]_B_X + mesh[i + 1][j]_B_Y * mesh[i + 1][j]_B_Y +
                         mesh[i + 1][j]_B_Z * mesh[i + 1][j]_B_Z;
            double b23 = mesh[i][j + 1]_B_X * mesh[i][j + 1]_B_X + mesh[i][j + 1]_B_Y * mesh[i][j + 1]_B_Y +
                         mesh[i][j + 1]_B_Z * mesh[i][j + 1]_B_Z;
            double b24 =
                    mesh[i + 1][j + 1]_B_X * mesh[i + 1][j + 1]_B_X + mesh[i + 1][j + 1]_B_Y * mesh[i + 1][j + 1]_B_Y +
                    mesh[i + 1][j + 1]_B_Z * mesh[i + 1][j + 1]_B_Z;

            double gammag = PhysConsts::gamma;
            double c1 = (gammag * pressurearr[i][j] + b21) * rho1i;
            double c2 = (gammag * pressurearr[i + 1][j] + b22) * rho2i;
            double c3 = (gammag * pressurearr[i][j + 1] + b23) * rho3i;
            double c4 = (gammag * pressurearr[i + 1][j + 1] + b24) * rho4i;

            double v1 = mesh[i][j]_MOMX / mesh[i][j]_MASS;
            double v2 = mesh[i][j + 1]_MOMX / mesh[i][j + 1]_MASS;
            double v3 = mesh[i + 1][j]_MOMX / mesh[i + 1][j]_MASS;
            double v4 = mesh[i + 1][j + 1]_MOMX / mesh[i + 1][j + 1]_MASS;

            double v5 = mesh[i + 1][j - 1]_MOMY / mesh[i + 1][j - 1]_MASS;
            double v6 = mesh[i][j + 1]_MOMY / mesh[i][j + 1]_MASS;
            double v7 = mesh[i + 1][j]_MOMY / mesh[i + 1][j]_MASS;
            double v8 = mesh[i + 1][j + 1]_MOMY / mesh[i + 1][j + 1]_MASS;

            double dv1 = (mesh[i + 1][j + 1]_MOMX - mesh[i][j + 1]_MOMX);
            double dv2 = (mesh[i + 1][j]_MOMX - mesh[i][j]_MOMX);
            double dv3 = (mesh[i + 1][j + 1]_MOMY - mesh[i + 1][j]_MOMY);
            double dv4 = (mesh[i + 1][j]_MOMY - mesh[i + 1][j - 1]_MOMY);

            double divu = 0.5 * (dv1 - dv2 + dv3 - dv4);
            if (-betaB * min4(c1, c2, c3, c4) > delta_x * (divu)) {
                swb = true;
            }

            double psi = 0;
            /*
            if ( swa && swb)
            {
                psi=delpz /(delpz + delpr);
                std::cout << "Upwinding Ez" << std::endl;
            }
            else {
                psi=0.5;
            }
            */
            psi = 0.5;

            Ez[i][j] = 0.0;
            Ez[i][j] =
                    0.5 * psi * (-faceEx[i][j] - faceEx[i][j - 1])
                    + 0.5 * (1.0 - psi) * (faceEy[i][j] + faceEy[i - 1][j]);
            //Ez[i][j] =1.0;

//      double oldway = 0.25 * (-xflux[i][j] _B_Y - xflux[i][j - 1] _B_Y + yflux[i][j] _B_X + yflux[i - 1][j] _B_X);
//      Ez[i][j] = oldway;
#ifdef OLDWAY_CHECK22
            if (oldway != Ez[i][j])
              {
                std::cout
                  << "Error: fluxes not equal face fields"
                  << " ( " << i << " , " << j << " ) " << std::endl;
                MPI_Abort (MPI_COMM_WORLD, 2222222);
                exit (0);
              }
#endif

#ifdef  RESISTIVE



            double myx = (double) myaddress[0] + i -NGC+0.5;
            double myy = (double) myaddress[1] + j -NGC +0.5;
            double r = myx * delta_x;
            double z = myy * delta_y;

            double eps = 0.1;
            double GM = 1.0 / (eps * eps);
            double etaz = 0.1;
            double etam = 0.1;
            double alpham = 0.1;
            // H is the thermal heightscale of the disk
            // cs is soundspeed at disk midplane
            // cs = sqrt (gamma *P/ rho);
            double cs = 0;
            cs = SoundSpeedMidplane[i];
            // omegaK is keplerian velocity at disk midplane
            double omegaK = 0;
            omegaK = sqrt (GM / r * r * r);
            double H = eps * r;
            H = cs / omegaK;
            double bx = mesh[i][j] _B_X;
            double by = mesh[i][j] _B_Y;
            double bz = mesh[i][j] _B_Z;
            double Va =
              std::sqrt (bx * bx + by * by +
                     bz * bz) / std::sqrt (mesh[i][j] _MASS);
            etam = alpham * H * Va * exp (-2 * z * z / (H * H));
            double chim = 1.0;
            double etamdash = chim * etam;
            double etaR = etamdash;
            double etaZ = etamdash;
            double etaPhi = etam;
            Ez[i][j] +=
              0.25 * etaZ * (xjz[i][j] - xjz[i - 1][j] + yjz[i][j] - yjz[i][j - 1]);
              double aa=Ez[i][j];
              if ( fabs(aa)> 1e10)
              {
                  std::cout
                  << "Nonphysical Ez = "  << aa
                  << " ( " << i << " , " << j << " ) "
                  << " " << etaZ
                  << " alpham " << alpham
                  << " H " << H
                  << " Va " << Va
                  << " bx " << bx
                  << " by " << by
                  << " bz " << bz
                   << " " << xjz[i][j] << " " <<  xjz[i - 1][j] << " " <<  yjz[i][j] << " " <<  yjz[i][j - 1]
                  << std::endl;
                  exit(0);
              }
#endif
        }


#ifdef  STAGGER_MESH

/*
  for (int i = 2; i < nx - 1; ++i)
  {
    for (int j = 2; j < ny - 1; ++j)
		{

			double tst= Ez[i][j];
			if (fabs (tst) > 1e-14 ) { 
			std::cout << tst << std::endl;
			std::cout 
			<< i  << " "  << j 
			<< std::endl;

			std::cout 
	      << faceEx[i][j] 
			<< " " 
			<< faceEx[i][j - 1] 
			<< " " 
			<< faceEy[i][j] 
			<< " " 
			<< faceEy[i - 1][j]
			<< std::endl;
			exit(0);
			}
		}
  }
  */


/*
  for (int i = 2; i < ny ; ++i)
    {
//      Ez[2][j] = 0.5 * ( faceEy[2][j] + faceEy[2 - 1][j]);
      Ez[i][0] =   Ez[i][2];
      Ez[i][1] =   Ez[i][3];
//      Ez[i][ny-2] =  -Ez[i][ny-3];
//      Ez[i][ny-1] =  -Ez[i][ny-2];
    }

  for (int j = 0; j < ny ; ++j)
    {
//      Ez[2][j] = 0.5 * ( faceEy[2][j] + faceEy[2 - 1][j]);
      Ez[0][j] =  Ez[4][j];
      Ez[1][j] =  Ez[3][j];
      Ez[2][j] = 0.;
//      Ez[ny-2][j] = Ez[ny-3][j]  ;
 //     Ez[ny-1][j] = Ez[ny-2][j]  ;
    }
	 */

    // Advect face-centred fields using Faraday's law
    for (int i = 2; i < nx - 1; ++i)
        for (int j = 2; j < ny - 1; ++j) {
            faceBx[i][j] -= (dtody * (Ez[i][j + 1] - Ez[i][j]));

            double myx = (double) myaddress[0] + (i - NGC);
            double dri = 2.0 / (2.0 * myx + 1);
#ifdef CYLINDRICAL
            // Cylindrical curl operator NRL Formulary page 7
            //faceBy[i][j + 1] += dtodx * ((myx + 1.0) * Ez[i + 1][j + 1] - (myx - 0.0) * Ez[i][j + 1]) * dri;
            faceBy[i][j ] += dtodx * ((myx + 1.0) * Ez[i + 1][j ] - (myx - 0.0) * Ez[i][j]) * dri;
            //faceBy[i][j ] += dtodx * ((myx + 0.5) * Ez[i + 1][j ] - (myx - 0.5) * Ez[i][j]) * dri;
#else
            faceBy[i][j] += (+dtodx) * (Ez[i + 1][j] - Ez[i][j]);
#endif

#undef DEBUGG
#ifdef DEBUGG
            std::cout
            <<  "xyxy"
            <<  " " << i
            <<  " " << j
            <<  " " << faceBx[i][j ]
            <<  " " << Ez[i][j+1 ]
            <<  " " << Ez[i][j ]
            << std::endl
            <<  " " << faceBy[i][j ]
            <<  " " << Ez[i+1][j ]
            <<  " " << Ez[i][j ]
            << std::endl
            ;
#endif
        }
#endif

    std::cout << faceBx[5][5] << std::endl;


#ifdef  STAGGER_MESH
#define CORRECTION
#endif

#ifndef CORRECTION
    if (myid == 1) { std::cout << "No field correction" << std::endl; }
#endif

#ifdef CORRECTION
    for (int i = 2; i < nx - 2; ++i)
        for (int j = 2; j < ny - 2; ++j) {
            // Correct the Energy to account for the staggered mesh
            double et = newmesh[i][j]_ENER;
            double bx = newmesh[i][j]_B_X;
            double by = newmesh[i][j]_B_Y;
            double bz = newmesh[i][j]_B_Z;
            double centredb2 = 0.5 * (bx * bx + by * by);
            double avebx = 0.5 * (faceBx[i + 1][j] + faceBx[i][j]);

#ifdef CYLINDRICAL
            //double myx = myaddress[0] + (i - NGC);
            //double rgi = (myx * myx + myx + 1.0 / 3.0) * delta_x / (myx + 0.5);
            //rgi = (myx * myx - myx + 1.0 / 3.0) * delta_x / (myx - 0.5);
           // avebx= faceBx[i + 1][j] - idx * ((i + 1) * delta_x - rgi) * (faceBx[i + 1][j] - faceBx[i][j]);
              avebx = 0.5 * (faceBx[i + 1][j] + faceBx[i][j]);
#endif
            double aveby = 0.5 * (faceBy[i][j + 1] + faceBy[i][j]);
            double faceb2 = 0.5 * (avebx * avebx + aveby * aveby);
            double correction = faceb2 - centredb2;
            if (correction > 0) {
                newmesh[i][j]_ENER += correction;
            }


#ifdef CHANGE_E
            double rho, rhoi;
            double px;
            double py, pz, velx, vely, velz, velx2, vely2, velz2, ke, pressure;
            double gammag = PhysConsts::gamma;
            double gammam1 = gammag - 1;
            double gammam1i = 1.0 / gammam1;
            rho = newmesh[i][j] _MASS;
            px = newmesh[i][j].array[1];
            py = newmesh[i][j].array[2];
            pz = newmesh[i][j].array[3];
            et = newmesh[i][j].array[4];
            bx = avebx;
            by = avebx;
            bz = mesh[i][j].array[7];
            rhoi = 1.0 / rho;
            velx = px * rhoi;
            vely = py * rhoi;
            velz = pz * rhoi;
            velx2 = velx * velx;
            vely2 = vely * vely;
            velz2 = velz * velz;
            ke = 0.5 * rho * (velx2 + vely2 + velz2);
            double bsquared = 0.5 * (bx * bx + by * by + bz * bz);
            double internalenergy = et - ke - bsquared;
            pressure = internalenergy * (gammam1);
#ifdef DONT
            if (pressure < 0)
              {
                std::cout << "Adj press" << std::endl;
                pressure = 0.25 * (std::max (pressurearr[i][j + 1], 0) +
                           std::max (pressurearr[i][j - 1], 0) +
                           std::max (pressurearr[i + 1][j], 0) +
                           std::max (pressurearr[i - 1][j], 0));
                pressurearr[i][j] = pressure;
              }
#endif
            //pressure = std::max (pressure, pmin);
            et = ke + bsquared + pressure * gammam1i;
            newmesh[i][j] _ENER = et;
            //newmesh[i][j] _B_X = avebx;
            //newmesh[i][j] _B_Y = aveby;
#endif


        }
#endif

#ifdef GRAVITY
    for (int i = 2; i < nx - 2; i++)
      {
        for (int j = 2; j < ny - 2; j++)
      {
        // Gravity only in source term
        // Colella & Woodward 1984

        // Find distance
        //


        double myx = (double) myaddress[0] + i - NGC;
        double myy = (double) myaddress[1] + j - NGC;
        double r = myx * delta_x;
        double z = myy * delta_y;

        // rho g
        // rho u g
        double rho = mesh[i][j] _MASS;
        double eps = 0.1;
        double GM = 1.0 / (eps * eps);

        double g =
          (-GM) / std::sqrt (std::pow (r, 2.0) + std::pow (z, 2.0));
        newmesh[i][j] _MOMX += delta_t * rho * g;
        newmesh[i][j] _MOMY += delta_t * rho * g;
        newmesh[i][j] _ENER +=
          delta_t * (mesh[i][j] _MOMX * g + mesh[i][j] _MOMY * g);


      }
      }

#endif

// Check div B
#define CHECK_DIVB_HERE
#ifdef CHECK_DIVB_HERE
    TNT::Array2D<double> divb(nx, ny);
    divb = 0.0;
    double maxdivb = (-9.0e99);
    double mindivb = 9.0e99;
    for (int i = 2; i < nx - 2; i++) {
        for (int j = 2; j < ny - 2; j++) {

            double a = faceBx[i + 1][j];
            double b = faceBx[i][j];
            double c = faceBy[i][j + 1];
            double d = faceBy[i][j];

            divb[i][j] =
                    idx * (a - b) +
                    idy * (c - d);
//	  divb[i][j] = (faceBx[i + 1][j] - faceBx[i][j]) + (faceBy[i][j + 1] - faceBy[i][j]);
#ifdef CYLINDRICAL
            double myx = (double) myaddress[0] + (i - NGC);
            double dri = 2.0 / (2 * myx + 1);
            divb[i][j] =
              idx * dri * ((myx + 1.0) * faceBx[i + 1][j] -
                   (myx - 0.0) * faceBx[i][j]) + idy * (faceBy[i][j +
                                          1] -
                                        faceBy[i][j]);
#endif
            newmesh[i][j]_DIVB = divb[i][j];

            maxdivb = std::max(maxdivb, divb[i][j]);
            mindivb = std::min(mindivb, divb[i][j]);
            if (0 && fabs(divb[i][j]) > 1e-20)
                //if (fabs (divb[i][j]) > 1e-20)
            {

                std::cout << std::setprecision(25);
                std::cout
                        << " old " <<
                        idx * (oldfaceBx[i + 1][j] - oldfaceBx[i][j]) +
                        idy * (oldfaceBy[i][j + 1] - oldfaceBy[i][j])
                        << std::endl;


                std::cout
                        << "(" << i << "," << j << ")" << divb[i][j] << std::endl;

                std::cout
                        << " " << faceBx[i + 1][j]
                        << " " << oldfaceBx[i + 1][j]
                        << " " << dtody
                        << " " << Ez[i + 1][j + 1] << " " << Ez[i + 1][j]
                        << std::endl;
                std::cout
                        << " " << faceBx[i][j]
                        << " " << oldfaceBx[i][j]
                        << " " << dtody
                        << " " << Ez[i][j + 1] << " " << Ez[i][j] << std::endl;
                std::cout
                        << " " << faceBy[i][j + 1]
                        << " " << oldfaceBy[i][j + 1]
                        << " " << dtodx
                        << " " << Ez[i + 1][j + 1] << " " << Ez[i][j + 1]
                        << std::endl;
                std::cout
                        << " " << faceBy[i][j]
                        << " " << oldfaceBy[i][j]
                        << " " << dtodx
                        << " " << Ez[i + 1][j] << " " << Ez[i][j] << std::endl;
                std::cout
                        << -dtody * (Ez[i + 1][j + 1] - Ez[i + 1][j]) +
                           dtody * (Ez[i][j + 1] - Ez[i][j]) +
                           dtodx * (Ez[i + 1][j + 1] - Ez[i][j + 1]) -
                           dtodx * (Ez[i + 1][j] - Ez[i][j]) << std::endl;
                std::cout
                        << "Bombing out for div b"
                        << std::endl;
                exit(0);
            }
        }
    }

    std::cout << std::setiosflags(std::ios::scientific);
    std::cout << " MaxDivB = " << maxdivb << std::endl;
    std::cout << " MinDivB = " << mindivb << std::endl;
#endif


    for (int i = 2; i < nx - 2; i++) {
        for (int j = 2; j < ny - 2; j++) {
            //final pressure check


            double et = newmesh[i][j]_ENER;
            double bx = newmesh[i][j]_B_X;
            double by = newmesh[i][j]_B_Y;
            bx = 0.5 * (faceBx[i][j] + faceBx[i + 1][j]);
            by = 0.5 * (faceBy[i][j] + faceBy[i][j + 1]);
            double bz = newmesh[i][j]_B_Z;
            double rho, rhoi;
            double px;
            double py, pz, velx, vely, velz, velx2, vely2, velz2, ke, pressure;
            double gammag = PhysConsts::gamma;
            double gammam1 = gammag - 1;
            double gammam1i = 1.0 / gammam1;
            rho = newmesh[i][j]_MASS;
            px = newmesh[i][j].array[1];
            py = newmesh[i][j].array[2];
            pz = newmesh[i][j].array[3];
            bz = newmesh[i][j].array[7];
            rhoi = 1.0 / rho;
            velx = px * rhoi;
            vely = py * rhoi;
            velz = pz * rhoi;
            velx2 = velx * velx;
            vely2 = vely * vely;
            velz2 = velz * velz;
            ke = 0.5 * rho * (velx2 + vely2 + velz2);
            double bsquared = 0.5 * (bx * bx + by * by + bz * bz);
            double internalenergy = et - ke - bsquared;
            pressure = internalenergy * (gammam1);
            if (pressure < 0) {
                //	newmesh[i][j] _ENER =  ke + bsquared + pmin*gammam1i;
#define SADSTORY
#ifdef SADSTORY
                std::cout << "Neg pressure = "
                          << pressure << " "
                          << " at (" << i << "," << j << ")" << std::endl;

                double g1 = 0.5 * (mesh[i][j]_MOMX * mesh[i][j]_MOMX + mesh[i][j]_MOMY * mesh[i][j]_MOMY +
                                   mesh[i][j]_MOMZ * mesh[i][j]_MOMZ) /
                            mesh[i][j]_MASS;
                double g2 = 0.5 * (mesh[i][j]_B_X * mesh[i][j]_B_X + mesh[i][j]_B_Y * mesh[i][j]_B_Y +
                                   mesh[i][j]_B_Z * mesh[i][j]_B_Z);

                std::cout
                        << " " << (mesh[i][j]_ENER - g1 - g2) * gammam1
                        << " " << mesh[i][j]_ENER
                        << " " << g1
                        << " " << g2
                        << std::endl;
                std::cout
                        << pressure
                        << " = " << et
                        << " - " << ke
                        << " - " << bsquared
                        << std::endl;
                // What is the sad story behind these negative pressures?
                for (int hh = 0; hh < 8; ++hh) {
                    std::cout << " " << mesh[i + 1][j].array[hh];
                }
                std::cout << std::endl;
                for (int hh = 0; hh < 8; ++hh) {
                    std::cout << " " << mesh[i - 1][j].array[hh];
                }
                std::cout << std::endl;
                for (int hh = 0; hh < 8; ++hh) {
                    std::cout << " " << mesh[i][j + 1].array[hh];
                }
                std::cout << std::endl;
                for (int hh = 0; hh < 8; ++hh) {
                    std::cout << " " << mesh[i][j - 1].array[hh];
                }
                std::cout << std::endl;
                std::cout
                        << " newmesh"
                        << "\tmesh"
                        << "\tdtodx"
                        << std::endl;
                for (int hh = 0; hh < 8; ++hh) {
                    std::cout
                            << " " << newmesh[i][j].array[hh]
                            << " " << mesh[i][j].array[hh]
                            << " " << dtodx
                            << " " << xflux[i + 1][j].array[hh]
                            << " " << xflux[i][j].array[hh]
                            << " " << dtody
                            << " " << yflux[i][j + 1].array[hh]
                            << " " << yflux[i][j].array[hh] << std::endl;
                }
                std::cout
                        << " " << bx
                        << "  " << by
                        << std::endl;
#endif

            }

        }
    }

/*    int q1 = 202, q2 = 2;
    std::cout << "faceBx"
              << " " << q1
              << " " << q2
              << " " << faceBx[q1][q2]
              << " " << Ez[q1][q2 + 1]
              << " " << Ez[q1][q2] << std::endl;

    std::cout << "Ez "
              << " " << Ez[q1][q2 + 1]
              << " " << faceEx[q1][q2 + 1]
              << " " << faceEx[q1][q2]
              << " " << faceEy[q1][q2 + 1]
              << " " << faceEy[q1 - 1][q2 + 1]
              << std::endl;*/

#ifdef COOLING
    for (int i = 2; i < nx - 2; i++)
      {
        for (int j = 2; j < ny - 2; j++)
      {

          double Lcool=0;
          cooling (mesh[i][j].array, &Lcool,delta_t);
          mesh[i][j] _ENER +=Lcool;
          mesh[i][j] _COOL =Lcool;
      }
       }
#endif
    return 0;
}
