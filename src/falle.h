/* $Id: falle.h,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "global.h"
#include "eigenvectors.h"
int falle_mhd (double *leftstate,
	       double *rightstate,
	       double *, int time_step, double *max_speed, int idir);
