/* $Id: mhd_godunov.h,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
/* $Id: mhd_godunov.h,v 1.3 2006-10-30 15:17:55 gmurphy Exp $: mhd_godunov.h,v 1.2 2006-10-30 15:12:38 gmurphy Exp $ */
//#include <math.h>
//#include <stdio.h>
//#include "burger.h"
#include "global.h"
#include "falle.h"
int godunov_mhd (zone ** grid, zone ** flux1, zone ** flux2,
		 int height, int width, int jj, double *time);
