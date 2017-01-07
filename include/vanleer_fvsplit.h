/* $Id: vanleer_fvsplit.h,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "global.h"
int vanleer_flux_vector_split (double *leftstate, double *rightstate,
			       double *fp, double *fn, int iii, int jjj,
			       int idir);
