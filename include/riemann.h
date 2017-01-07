/* $Id: riemann.h,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "global.h"
#include "eigenvectors.h"
int riemann (double *leftstate,
	     double *rightstate,
	     double *roeflux,
	     double *res_state, int time_step, double *max_speed, int idir);
