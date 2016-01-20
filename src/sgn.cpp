/* $Id: sgn.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "sgn.h"

double
sgn (double a)
{
  if (a < 0)
    {
      return -1.0;
    }
  else
    {
      return 1.0;
    }
}
