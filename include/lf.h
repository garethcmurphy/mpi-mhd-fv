/* $Id: hlld.h,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "global.h"
#include "eigenvectors.h"
int lf (double *leftstate,
	  double *rightstate,
	  double *fhlld,
	  double *Res_state,
	  int time_step, double *max_speed, int idir, int *loc);
