/* $Id: maes.cpp,v 1.5 2006-11-16 13:48:07 gmurphy Exp $  */


#include "initialise_maes.h"

int
zannisimple(TNT::Array2D<unk> mesh,
            TNT::Array2D<double> faceBx,
            TNT::Array2D<double> faceBy,
            MPI_Comm new_comm, int ndims, int *dim_size,
            double *xsize, double *ysize) {
    double minpressure = 9e99;
    double mindensity = 9e99;
    double maxdensity = (-9e99);

    std::cout << "Initialise Magnetised Accretion-Ejection Structure " << std::
    endl;

    const int nx = mesh.dim1();
    const int ny = mesh.dim2();
    const double gammam1i = 1 / (PhysConsts::gamma - 1);
    const double gammam1 = (PhysConsts::gamma - 1);
    double r, z, h;
    double r2, z2, h2;

// Normalization
// r0 = 0.1 AU for YSO
// r0 = 10 RSchw for AGN
    const double r0 = 1.0;
    // Disk Aspect Ratio
    const double eps = 0.1;
    const double vk0 = 1.0;
    const double rho0 = 1.0;
    const double cs0 = eps * vk0;
    const double pres0 = rho0 * cs0 * cs0;
    // The Zanni et al 6 free parameters

    // Disk magnetization parameter
    const double mu = 0.3;
    // Magnetic heightscale
    // The parameter m determines the heightscale on which the initial magnetic
    // field bends
    // m -> inf is perfectly vertical
    //
    const double m = 0.35;

    double r02;
    // ms is a parameter smaller than unity
    // which then ensures an initial subsonic poloidal
    // inflow
    double ms = 0.3;
    double rho, vtheta, vz, vr, bz, br, btheta;
    // plasma beta parameter measuring the ratio of the
    // thermal pressure to the magnetic pressure at Z = 0
    double beta = 1.0;

    double pressure = 0;

    int coords[] = {0, 0};

    int myid = 99;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    MPI_Cart_coords(new_comm, myid, ndims, coords);
    int dims[] = {0, 0};
    MPI_Cartdim_get(new_comm, dims);

    int xmax = 0;
    int ymax = 0;
    MPI_Allreduce(coords, &xmax, 1, MPI_INT, MPI_MAX, new_comm);
    MPI_Allreduce(coords + 1, &ymax, 1, MPI_INT, MPI_MAX, new_comm);
    xmax++;
    ymax++;

// Get my address in the global sim
    std::cout << "Proc" << myid << " " << xmax << " " << ymax << std::endl;
    int myaddress[] = {0, 0};
    myaddress[0] = (nx - 2 * NGC) * coords[0];
    myaddress[1] = (ny - 2 * NGC) * coords[1];


    double bz0 = 1.0;

// Set the size of the simulation
    *xsize = 20 * r0;
    *ysize = 20 * r0;

    double delta_x = *xsize / ((nx - 2 * NGC) * xmax);
    double delta_y = *ysize / ((ny - 2 * NGC) * ymax);

    TNT::Array2D<double> Psi(nx + 1, ny + 1);

    // Initialise the corner-centred magnetic flux function
    // according to Zanni et al 2007, Casse-Keppens 2004
    for (int jj = 0; jj < ny + 1; jj++)
        for (int ii = 0; ii < nx + 1; ii++) {
            double myx = 0;
            double myy = 0;
            myx = (double) myaddress[0] + ii - NGC;
            myy = (double) myaddress[0] + jj - NGC;

            r = (*xsize * myx) / (dim_size[0] * (nx - 2 * NGC)) + 1e-20;
            z = (*ysize * myy) / (dim_size[1] * (ny - 2 * NGC));
            z2 = z * z;
            r2 = r * r;

/*
	if (r == 0)
	{
	Psi[ii][jj] = 0;
	}
	else
	*/
            {
                Psi[ii][jj] = (4.0 / 3.0) * (bz0 * r0 * r0) * pow(r / r0, 0.75) * pow(r, 1.25) * pow(m, 1.25) /
                              pow((m * m * r2 + z2), 5.0 / 8.0);
                Psi[ii][jj] = (4.0 / 3.0) * (bz0 * r0 * r0) * pow(r / r0, 0.75) * pow(m, 1.25) /
                              pow((m * m + z2 / r2), 5.0 / 8.0);
                //Psi[ii][jj] =  r2;
            }

            //std::cout << Psi[ii][jj] << std::endl;

        }
    //


    // Initialise the face-centred magnetic field
    for (int jj = 2; jj < ny; jj++)
        for (int ii = 2; ii < nx; ii++) {
            //
            double myx = 0;
            double myy = 0;
            myx = (double) myaddress[0] + ii - NGC;
            //std::cout << delta_x << " " <<delta_y << std::endl;
            r = (*xsize * myx) / (dim_size[0] * (nx - 2 * NGC)) + delta_x * 0.5;
            h = eps * r;
            h2 = h * h;
            r2 = r * r;
            r02 = r0 * r0;
            bz = pow(r0, 2.5) / (pow((r02 + r2), 1.25) * sqrt(beta));
            faceBx[ii][jj] = 0.0;
            faceBy[ii][jj] = bz;

            double ri = 1 / r;
            double idx = 1 / delta_x;
            double idy = 1 / delta_y;
            double rxface = (*xsize * myx) / (dim_size[0] * (nx - 2 * NGC)) + 1e-20;
            double irxface = (1 / rxface);
            faceBx[ii][jj] = -irxface * idy * (Psi[ii][jj + 1] - Psi[ii][jj]);


/*
	std::cout 
	<< faceBx[ii][jj] 
	 << " " << -irxface
	 << " " << Psi[ii][jj+1]
	 << " " << Psi[ii][jj]
	 << std::endl;
	 */


            faceBy[ii][jj] = ri * idx * (Psi[ii + 1][jj] - Psi[ii][jj]);

/*
	std::cout 
	<< faceBy[ii][jj] 
	 << " dsize " << dim_size[0]
	 << " rxface " << rxface
	 << " " << Psi[ii+1][jj]
	 << " " << Psi[ii][jj]
	 << " " << rxface*rxface
	 << std::endl;
	 */

            if (isnan(bz)) {
                MPI_Abort(MPI_COMM_WORLD, 33333333);
                exit(0);
            }
//      faceBy[ii][jj] = 3.0;
        }

    for (int jj = 0; jj < ny + 1; jj++) {
        faceBy[nx][jj] =
        faceBy[nx - 1][jj] =
                faceBy[nx - 2][jj];

    }

#undef PROBS_WITH_B
#ifdef PROBS_WITH_B
    for (int hh = 0; hh < ny + 1; hh++)
      for (int gg = 0; gg < nx + 1; gg++)
        {
      faceBy[hh][gg] = 3.0;
        }
#endif



// Initial midplane density
//

    TNT::Array1D<double> density0(nx);
    TNT::Array1D<double> pressure0(nx);
    TNT::Array1D<double> temperature0(nx);
    TNT::Array2D<double> pressurefunc(nx, ny);
    for (int ii = 2; ii < nx; ii++) {
        int gg = ii - NGC;
        double myx = (double) myaddress[0] + gg;
        r = (*xsize * myx) / (dim_size[0] * (nx - 2 * NGC)) + delta_x * 0.5;
        r = std::max(r, r0);
        //double temp = 1.0e4*(eps/0.1)*(eps/0.1)*pow(10*r/PhysConsts::AU , -1);
        double rho = rho0 * pow(r / r0, -1.5);
        double pres = pres0 * pow(r / r0, -2.5);
        mesh[ii][2]_MASS = rho;
        density0[ii] = rho;
        pressure0[ii] = pres;
        pressurefunc[ii][2] = pres;
        mesh[ii][2]_MASS = density0[ii];
        //std::cout << " " << rho << " " << pres << std::endl;
    }


// Cal Runge Kutta function and work out sound speed, then
// use that to  calculate density and pressure

    /*
    for (int ii = 2; ii < nx; ii++)
  for (int jj = 3; jj < ny; jj++)
      {
            int gg = ii - NGC;
            int hh = jj - NGC;
            double myx = 0;
            double myy = 0;
            myx = (double) myaddress[0] + gg;
            myy = (double) myaddress[1] + hh;
            r = ( *xsize * myx) / (dim_size[0] * (nx-2*NGC) ) + delta_x * 0.5;
            z = ( *ysize * myy) / (dim_size[1] * (ny-2*NGC)) + delta_y * 0.5;
            z2 = z * z;
            r2 = r * r;
            double sphere_r= std::sqrt(r2+z2);

            if (r > r0)
            {

                if ( mesh[ii][jj-1] _MASS > 1e-4 )
                {
                double fc2 =  mesh[ii][jj-1] _MASS/rho0 * pow (r/r0, 3.0/2.0 ) ;
                double ggam= PhysConsts::gamma* (PhysConsts::gamma-1);
                fc2 = pow (fc2, ggam);
            fc2 =  rungekutta ( fc2 , z/r,  delta_y/r);
                if ( fc2 < 0 )
                {
                double rhoa0=1e-4;
                double rho = rhoa0*pow(r0/sphere_r, gammam1i) ;
            double gamp= PhysConsts::gamma/gammam1;
            double gampi= 1/gamp;
            double GM = 1/(eps*eps);
            double pressure=rhoa0 * gampi * GM * pow(r0/sphere_r , gamp)/ r0;
             mesh[ii][jj] _MASS= rho;
             pressurefunc[ii][jj] = pressure;
                }
                else{
                    // in the disk
            double density = rho0* pow(r/r0,-1.5)* pow(fc2, PhysConsts::gamma);
            double pressure = pres0* pow(r/r0,-2.5)* pow(fc2, PhysConsts::gamma-1);
             mesh[ii][jj] _MASS= density;
             pressurefunc[ii][jj] = pressure;
                }

                }
                else
                {
                double rhoa0=1e-4;
                double rho = rhoa0*pow(r0/sphere_r, gammam1i) ;
            double gamp= PhysConsts::gamma/gammam1;
            double gampi= 1/gamp;
            double GM = 1/(eps*eps);
            pressure=rhoa0 * gampi * GM * pow(r0/sphere_r , gamp)/ r0;
             mesh[ii][jj] _MASS= rho;
             pressurefunc[ii][jj] = pressure;
                }
            }
            else{
                double rhoa0=1e-4;
                double rho = rhoa0*pow(r0/sphere_r, gammam1i) ;
            double gamp= PhysConsts::gamma/gammam1;
            double gampi= 1/gamp;
            double GM = 1/(eps*eps);
            pressure=rhoa0 * gampi * GM * pow(r0/sphere_r , gamp)/ r0;
             pressurefunc[ii][jj] = pressure;
             mesh[ii][jj] _MASS= rho;
            }

#ifdef WRITEOUT
                std::cout
                << " " <<  ii
                << " " <<  jj
                << " " <<   mesh[ii][jj] _MASS << std::endl;
#endif


        }
        */


    for (int ii = 2; ii < nx; ii++)
        for (int jj = 2; jj < ny; jj++) {
            // Mass Density
            int gg = ii - NGC;
            int hh = jj - NGC;
            double myx = 0;
            double myy = 0;
            myx = (double) myaddress[0] + gg;
            myy = (double) myaddress[1] + hh;
            //std::cout << delta_x << " " <<delta_y << std::endl;
            z = (*ysize * myy) / (dim_size[1] * (ny - 2 * NGC)) + delta_y * 0.5;
            r = (*xsize * myx) / (dim_size[0] * (nx - 2 * NGC)) + delta_x * 0.5;
            r = std::max(r, delta_x * 0.5);
            z = std::max(z, delta_y * 0.5);
            h = eps * r;
            h2 = h * h;
            z2 = z * z;
            r2 = r * r;
            double sphere_r = std::sqrt(r2 + z2);
            /* R0 is a constant set equal to 4 in our runs.
             * This offset radius in the denominator makes the
             * density regular up to R = 0.
             */
            r02 = r0 * r0;
            rho =
                    std::max(1e-6, (pow(r0, 1.5) / pow((r0 * r0 + r * r), 0.75))) *
                    std::pow(std::max(1e-6, (1 - 0.5 * (gammam1) * z2 / h2)), gammam1i);
            //mesh[ii][jj] _MASS= (pow(r0,1.5)/pow((r0*r0+r*r),0.75))*pow (max(1e-6,  (1- 0.5* (gammam1) * z2/h2 )), gammam1i ) ;
            rho = std::max(1e-6, rho);
            /*
            std::cout
            << rho
            << " " << ii << " " << jj
            << " " << r
            << " " <<
              std::pow (std::max (1e-6, (1 - 0.5 * (gammam1) * z2 / h2)), gammam1i)

            << std::endl;
            */


            /*
               mesh[ii][jj] _MASS =  0.2 * (
               mesh[ii][jj] _MASS  +
               mesh[ std::min(ii+1,nx-1)    ][jj] _MASS  +
               mesh[ std::max(ii-1,0)       ][jj] _MASS  +
               mesh[ii                                                      ][std::min (jj+1, ny-1) ] _MASS  +
               mesh[ii                                                      ][std ::max(jj-1,0)] _MASS
               );
             */

            rho = mesh[ii][jj]_MASS;
            // Rotation Profile
            vtheta =
                    (1 -
                     eps * eps) * sqrt(r0) * exp(-2 * z2 / h2) / (eps *
                                                                  pow((r02 + r2),
                                                                      0.25));
            if (isnan(mesh[ii][jj]_MOMZ)) {
                std::cout << "z-momentum is nan" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 0);
                exit(0);
            }

            // Poloidal Velocity
            vr =
                    (-ms) * sqrt(r0) * exp(-2 * z2 / h2) / (pow((r02 + r2), 0.25));
            if (isnan(mesh[ii][jj]_MOMX)) {
                std::cout << "x-momentum is nan" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 99);
                exit(0);
            }

            vz = vr * z / r;
            if (isnan(mesh[ii][jj]_MOMY)) {
                std::cout << "y-momentum is nan" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 101);
                exit(0);
            }

            // Magnetic Field
            bz = 0.5 * (faceBy[ii][jj] + faceBy[ii][jj + 1]);
            br = 0.5 * (faceBx[ii][jj] + faceBx[ii + 1][jj]);
            btheta = 0;

            // Pressure = rho ^ gamma
            pressure = pow(rho, PhysConsts::gamma);
            pressure = pressurefunc[ii][jj];
            if (pressure < 0) {
                std::cout
                        << "Pressure < 0"
                        << std::endl;
                exit(0);
            }

            // Energy
            //mesh[ii][jj] _MASS = rho;
            mesh[ii][jj]_MOMZ = rho * vtheta;
            mesh[ii][jj]_MOMX = rho * vr;
            mesh[ii][jj]_MOMY = mesh[ii][jj]_MOMX * z / r;
            mesh[ii][jj]_B_X = br;
            mesh[ii][jj]_B_Y = bz;
            mesh[ii][jj]_B_Z = 0.;
            double vv2 = vz * vz + vr * vr + vtheta * vtheta;
            double b2 = bz * bz + br * br + btheta * btheta;
            mesh[ii][jj]_ENER = (0.5 * rho * vv2) + (0.5 * b2) + pressure * gammam1i;
            if (mesh[ii][jj]_ENER < 0) {

                std::cout
                        << " Neg init energy"
                        << " (" << ii << " , " << jj << " )" << std::endl;
            }

            mesh[ii][jj].temperature = myid;

            double rhoa0 = 1e-4;
            double atmosphererho = rhoa0 * pow(r0 / sphere_r, gammam1i);
            double gamp = PhysConsts::gamma / gammam1;
            double gampi = 1 / gamp;
            double GM = 1 / (eps * eps);
            double atmospherepressure = rhoa0 * gampi * GM * pow(r0 / sphere_r, gamp) / r0;
            //pressure=pressurefunc[ii][jj];

            // DISK STUFF
            // SK = sub keplerian
            double SK = 1 - eps * eps * gamp;
            double func = (1 / sphere_r - SK / r) * (gammam1) / (PhysConsts::gamma * eps * eps);
            func = std::max(func, 0.0);
            double diskpressure = pres0 * pow(func, gamp);
            double diskdensity = pow(func, gammam1i);
            vtheta = std::sqrt(SK / r);

            if (diskpressure > atmospherepressure) {
                rho = diskdensity;
                pressure = diskpressure;
            } else {
                rho = atmosphererho;
                pressure = atmospherepressure;
                vtheta = 0;

            }

            mesh[ii][jj]_MASS = rho;
            mesh[ii][jj]_MOMZ = rho * vtheta;
            mesh[ii][jj]_MOMX = rho * vr;
            mesh[ii][jj]_MOMY = mesh[ii][jj]_MOMX * z / r;
            mesh[ii][jj]_MOMX = 0;
            mesh[ii][jj]_MOMY = 0;
            vv2 = vtheta * vtheta;
            mesh[ii][jj]_ENER = (0.5 * rho * vv2) + (0.5 * b2) + pressure * gammam1i;


            minpressure = std::min(pressure, minpressure);
            mindensity = std::min(rho, mindensity);
            maxdensity = std::max(rho, maxdensity);
        }


    std::cout << "Min pressure " << minpressure << std::endl;
    std::cout << "Min density " << mindensity << std::endl;
    std::cout << "Max density " << maxdensity << std::endl;

    //

    return 0;


}
