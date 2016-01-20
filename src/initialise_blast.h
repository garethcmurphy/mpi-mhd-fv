/* $Id: initialise_blast.h,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "out.h"
int initialise_blast (TNT::Array2D < unk > mesh,
		      TNT::Array2D < double >faceBx,
		      TNT::Array2D < double >faceBy,
		      MPI_Comm new_comm, int ndims, int *dim_size,
		      double *xsize, double *ysize);
