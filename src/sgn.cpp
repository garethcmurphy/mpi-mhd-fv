#include "sgn.h"

double
sgn(double a) {
    if (a < 0) {
        return -1.0;
    } else {
        return 1.0;
    }
}
