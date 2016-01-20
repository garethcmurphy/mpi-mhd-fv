/* $Id: emf_exchange.h,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */

#include "problem.h"
int
emf_exchange (TNT::Array2D < double >emfx,
	      TNT::Array2D < double >emfy,
	      int myNorth, int mySouth, int myEast, int myWest, int myid);
