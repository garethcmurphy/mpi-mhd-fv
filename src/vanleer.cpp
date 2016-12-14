/* $Id: vanleer.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include "vanleer.h"

double
vanleer(double a, double b) {
    double temp = 0;
    temp = a * a + b * b;
    if (temp > 0.0)
        return (a * a * b + b * b * a) / temp;
    else
        return 0;
}
