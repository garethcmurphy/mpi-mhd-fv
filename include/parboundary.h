/* $Id: parboundary.h,v 1.4 2006-11-16 13:48:07 gmurphy Exp $  */

#include "out.h"

int parboundary (
TNT::Array2D < unk > mesh,
TNT::Array2D < unk > newmesh,
		 TNT::Array2D < double >faceBx,
		 TNT::Array2D < double >faceBy,
		 int myNorth, int mySouth, int myEast, int myWest, int myid,
		 MPI_Comm Cart_comm);
