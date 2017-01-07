/* $Id: ctop.h,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */

#include "mpi.h"
#include "out.h"
#include "PhysConsts.h"
#include "problem.h"

int ctop (double *c, double *p, double *loc);
int ptoc (double *p, double *c);
