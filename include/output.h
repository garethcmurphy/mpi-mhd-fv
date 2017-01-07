/* $Id: output.h,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "global.h"
int output (Array3D < zone > grid, Array3D < zone > fx, Array3D < zone > fy,
	    int time, char *filename);
