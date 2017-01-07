/* $Id: initialise_jet.h,v 1.1 2006-11-05 12:51:25 gmurphy Exp $  */
#include "out.h"
int initialise_jet (TNT::Array2D < unk > mesh,
		    TNT::Array2D < double >faceBx,
		    TNT::Array2D < double >faceBy,
		    MPI_Comm new_comm, int ndims, int *dim_size,
		    double *xsize, double *ysize);
