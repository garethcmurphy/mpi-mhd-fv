/* $Id: initialise.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "initialise_maes.h"

int
initialise (char *filename, TNT::Array3D < zone > grid, int *maxstep, double *cfl)
{
  int ii = 0;
  int jj = 0;
  int kk = 0;
  double rr, vx, vy, vz, p, ke, b1, b2, b3;
  double rl, ul, vl, wl, pl, bul, bvl, bwl, al;
  double rho, ur, vr, wr, pr, bur, bvr, bwr, ar;
  ifstream input;


  int tempint;
  double tempdouble;

  double gammam1 = gammag - 1;
  double gammam1i = 1.0 / gammam1;
  double dist;

  double bsquared = 0;
  double sqr4pie = 0;
  double sqr4piei = 0;
  double unused = 0;
  sqr4pie = 2.0 * sqrt (3.14159);
  sqr4piei = 1.0 / sqr4pie;


  rl = 1.0;
  ul = 0.0;
  vl = 0.0;
  pl = 6.0;

  rr = 0.1;
  ur = 0.0;
  vr = 0.0;
  pr = 0.6;

  input.open (filename);

  input >> tempint;
  input.ignore (256, '\n');
  input >> tempint;
  input.ignore (256, '\n');
  input >> tempint;
  input.ignore (256, '\n');
  input >> tempdouble;
  input.ignore (256, '\n');
  input >> tempint;
  input.ignore (256, '\n');
  input >> tempdouble;
  input.ignore (256, '\n');
  input >> tempdouble;
  input.ignore (256, '\n');
  input >> tempint;
  input.ignore (256, '\n');
  input.ignore (256, '\n');
  input >> rl;
  input.ignore (256, '\n');
  input >> ul;
  input.ignore (256, '\n');
  input >> vl;
  input.ignore (256, '\n');
  input >> wl;
  input.ignore (256, '\n');
  input >> pl;
  input.ignore (256, '\n');
  input >> bul;
  input.ignore (256, '\n');
  input >> bvl;
  input.ignore (256, '\n');
  input >> bwl;
  input.ignore (256, '\n');

  input >> rr;
  input.ignore (256, '\n');
  input >> ur;
  input.ignore (256, '\n');
  input >> vr;
  input.ignore (256, '\n');
  input >> wr;
  input.ignore (256, '\n');
  input >> pr;
  input.ignore (256, '\n');
  input >> bur;
  input.ignore (256, '\n');
  input >> bvr;
  input.ignore (256, '\n');
  input >> bwr;
  input.ignore (256, '\n');

  input.close ();


  al = sqrt (gammag * pl / rl);
  ar = sqrt (gammag * pr / rr);


  cout << "Left state rho= " << rl << "\t Rightstate " << rr << endl;
  cout << "Left state vx = " << ul << "\t Rightstate " << ur << endl;
  cout << "Left state vy = " << vl << "\t Rightstate " << vr << endl;
  cout << "Left state vz = " << wl << "\t Rightstate " << wr << endl;
  cout << "Left state p  = " << pl << "\t Rightstate " << pr << endl;
  cout << "Left state Bx = " << bul << "\t Rightstate " << bur << endl;
  cout << "Left state By = " << bvl << "\t Rightstate " << bvr << endl;
  cout << "Left state Bz = " << bwl << "\t Rightstate " << bwr << endl;
  cout << "Left state Al = " << al << "\t Rightstate " << ar << endl;
  cout << endl;


  /* Initialise the grid array with a shock tube problem  */
  for (ii = 0; ii < nx; ii++)
    {
      for (jj = 0; jj < ny; jj++)
	{
	  vx = ul;
	  vy = vl;
	  vz = wl;
//          p = 1+ 0.1*sin( 20*3.14*(double)jj /100);
	  p = pl;
	  rho = rl;
	  b1 = bul;
	  b2 = bvl;
	  b3 = bwl;

//   JET
//   if ( jj < ny * 0.25 && ii < nx * 0.75 && ii > nx * 0.25)
//   SHOCK TUBE
//     if ( ii < (14) )
	  if (ii > nx * 0.5)
//       if ( jj > ny*0.5)
//        if (ii <4 && jj<4)
//   SHOCK TUBE
//       if ((ii - jj) < 0)
	    //      if ((ii + jj) > nx)
//   SHOCK TUBE
//       if ((ii + 0.5*jj) < nx*0.75)
//        if (ii >=49 && ii <=50 && jj>=49 && jj<=50)
//        if (ii >=49 && ii <=50 )
//
	    //     dist = sqrt ((double)( ii*ii+jj*jj));
	    //            if (dist < 0.5*nx ) 

	    {
	      vx = ur;
	      vy = vr;
	      vz = wr;
	      p = pr;
	      rho = rr;
	      b1 = bur;
	      b2 = bvr;
	      b3 = bwr;
	    }

	  b1 = b1 * sqr4piei;
	  b2 = b2 * sqr4piei;
	  b3 = b3 * sqr4piei;

	  bsquared = (b1 * b1 + b2 * b2 + b3 * b3);
	  grid[ii][jj][kk] _MASS = rho;
	  grid[ii][jj][kk] _MOMX = rho * vx;
	  grid[ii][jj][kk] _MOMY = rho * vy;
	  grid[ii][jj][kk] _MOMZ = rho * vz;
	  ke = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
	  grid[ii][jj][kk] _ENER = ke + p * gammam1i + 0.5 * bsquared;
	  grid[ii][jj][kk] _B_X = b1;
	  grid[ii][jj][kk] _B_Y = b2;
	  grid[ii][jj][kk] _B_Z = b3;
	  if (grid[ii][jj][kk] _ENER < 0.0)
	    {
	      cout << "Wtf?" << endl;
	      exit (0);
	    }
	}

    }

  return 0;
}
