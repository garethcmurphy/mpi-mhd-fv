/* $Id: roe.h,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "global.h"
int roe (double *leftstate,
	 double *rightstate,
	 double *flux, double *fn, int iii, int jjj, int timestep, int idir);
