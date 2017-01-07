/* $Id: parupdate.h,v 1.4 2006-11-16 13:48:07 gmurphy Exp $  */

#include "out.h"
int parupdate (TNT::Array2D < unk > newmesh,
	       TNT::Array2D < unk > sourcemesh,
	       TNT::Array2D < unk > mesh,
	       TNT::Array2D < flux > fx,
	       TNT::Array2D < flux > fy,
	       TNT::Array2D < flux > SecondOrdCorr,
	       TNT::Array2D < double >faceBx,
	       TNT::Array2D < double >faceBy,
	       TNT::Array2D < double >xjx,
	       TNT::Array2D < double >xjy,
	       TNT::Array2D < double >xjz,
	       TNT::Array2D < double >yjx,
	       TNT::Array2D < double >yjy,
	       TNT::Array2D < double >yjz,
	       TNT::Array2D < double >faceEx,
	       TNT::Array2D < double >faceEy,
	       TNT::Array1D < double >SoundSpeedMidPlane,
	       double del,
	       double delta_x,
	       double delta_y, MPI_Comm Cart_comm, int second);