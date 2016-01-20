/* $Id: minmod.cpp,v 1.4 2006-10-30 15:17:55 gmurphy Exp $  */
/* $Id: minmod.cpp,v 1.4 2006-10-30 15:17:55 gmurphy Exp $: minmod.cpp,v 1.3 2006-10-30 15:12:38 gmurphy Exp $ */
#include "minmod.h"
#include <cmath>
double
minmod (double a, double b)
{
  /*int temp=copysignf (1,*a) + copysignf (1, *b);
     if ( fabs(*b) > fabs(*a) )
     return temp*fabs(*a)/2;
     else
     return temp*fabs(*b)/2; */
  if ((a * b) > 0)
    {
      double fa = std::fabs (a);
      double fb = std::fabs (b);

      return fa < fb ? a : b;
    }
  else
    return 0.0;
}
