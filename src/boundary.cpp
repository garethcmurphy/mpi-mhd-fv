/* $Id: boundary.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "boundary.h"
int
boundary (Array3D < zone > gridb)
{
  int ii = 0;
  int jj = 0;
  int kk = 0;
  int hh = 0;

//  cout << "bc:" << gridb[2][2][0]_MASS << endl;


  /* Boundary Conditions */
  /* at this point just copying solutions from
   * inside the grid to the outside - needs work */

  for (hh = 0; hh < ne; hh++)
    {
      gridb[1][jj][kk].array[hh] = gridb[2][jj][kk].array[hh];
      gridb[0][jj][kk].array[hh] = gridb[1][jj][kk].array[hh];
      // Upper y-boundary
      gridb[nx - 2][jj][kk].array[hh] = gridb[nx - 3][jj][kk].array[hh];
      gridb[nx - 1][jj][kk].array[hh] = gridb[nx - 2][jj][kk].array[hh];
    }

#ifdef TWODIM
  for (ii = 1; ii < nx - 1; ii++)
    {
      for (hh = 0; hh < ne; hh++)
	{
	  // Lower x-boundary
	  /* Copy row 2 to row 1 */
	  gridb[ii][1][kk].array[hh] = gridb[ii][2][kk].array[hh];
	  /* Copy row 1 to row 0 */
	  gridb[ii][0][kk].array[hh] = gridb[ii][1][kk].array[hh];
	  // Upper x-boundary
	  gridb[ii][ny - 2][kk].array[hh] = gridb[ii][ny - 3][kk].array[hh];
	  gridb[ii][ny - 1][kk].array[hh] = gridb[ii][ny - 2][kk].array[hh];
	}
    }

  for (jj = 1; jj < ny - 1; jj++)
    {
      for (hh = 0; hh < ne; hh++)
	{
	  // Lower y-boundary
	  gridb[1][jj][kk].array[hh] = gridb[2][jj][kk].array[hh];
	  gridb[0][jj][kk].array[hh] = gridb[1][jj][kk].array[hh];
	  // Upper y-boundary
	  gridb[nx - 2][jj][kk].array[hh] = gridb[nx - 3][jj][kk].array[hh];
	  gridb[nx - 1][jj][kk].array[hh] = gridb[nx - 2][jj][kk].array[hh];
	}
    }
#endif


#ifdef TWODIM
  // Corner values
  for (hh = 0; hh < ne; hh++)
    {
      gridb[1][1][kk].array[hh] = gridb[2][2][kk].array[hh];
      gridb[1][0][kk].array[hh] = gridb[2][0][kk].array[hh];
      gridb[0][1][kk].array[hh] = gridb[0][2][kk].array[hh];
      gridb[0][0][kk].array[hh] = gridb[1][1][kk].array[hh];

      gridb[nx - 2][1][kk].array[hh] = gridb[nx - 3][2][kk].array[hh];
      gridb[nx - 2][0][kk].array[hh] = gridb[nx - 3][0][kk].array[hh];
      gridb[nx - 1][1][kk].array[hh] = gridb[nx - 1][2][kk].array[hh];
      gridb[nx - 1][0][kk].array[hh] =
	0.5 * (gridb[nx - 2][0][kk].array[hh] +
	       gridb[nx - 1][1][kk].array[hh]);

      gridb[1][ny - 2][kk].array[hh] = gridb[2][ny - 3][kk].array[hh];
      gridb[1][ny - 1][kk].array[hh] = gridb[2][ny - 1][kk].array[hh];
      gridb[0][ny - 2][kk].array[hh] = gridb[0][ny - 3][kk].array[hh];
      gridb[0][ny - 1][kk].array[hh] =
	0.5 * (gridb[0][ny - 2][kk].array[hh] +
	       gridb[1][ny - 1][kk].array[hh]);
#ifdef DEBUG_BC
      cout << gridb[0][ny - 1][kk].array[hh] << endl;
#endif /*DEBUG*/
	gridb[nx - 2][ny - 2][kk].array[hh] =
	gridb[nx - 3][ny - 3][kk].array[hh];
      gridb[nx - 2][ny - 1][kk].array[hh] =
	gridb[nx - 3][ny - 1][kk].array[hh];
      gridb[nx - 1][ny - 2][kk].array[hh] =
	gridb[nx - 1][ny - 3][kk].array[hh];
      gridb[nx - 1][ny - 1][kk].array[hh] =
	gridb[nx - 2][ny - 2][kk].array[hh];

    }
#endif /*TWODIM*/
    /* End Boundary Conditions */
//  cout << "bc:" << gridb[2][2][0]_MASS << endl;
    return 0;
}
