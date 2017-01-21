//#include <math.h>
//#include <stdio.h>
//#include "burger.h"
#include "global.h"
#include "falle.h"
int godunov_mhd (zone ** grid, zone ** flux1, zone ** flux2,
		 int height, int width, int jj, double *time);
