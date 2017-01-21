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
