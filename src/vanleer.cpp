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
