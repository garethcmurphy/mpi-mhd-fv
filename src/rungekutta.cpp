

#include "PhysConsts.h"
#include <cmath>

double func(double x, double y);

double rungekutta(double init_y, double init_x, double h) {


    double k1 = h * func(init_x, init_y);
    double k2 = h * func(init_x + h / 2, init_y + k1 / 2);
    double k3 = h * func(init_x + h / 2, init_y + k2 / 2);
    double k4 = h * func(init_x + h, init_y + k3);


#ifdef XXX
    std::cout << init_x
    << " " << init_y
    << " " << init_y+(1.0/6.0)*(k1+k2+k3+k4)
    <<  std::endl;
#endif


    return init_y + (1.0 / 6.0) * (k1 + k2 + k3 + k4);
}

double func(double x, double y) {

    const double eps = 0.1;
    const double gameps = (1 - PhysConsts::gamma) / (PhysConsts::gamma * eps * eps);
    double tmp = gameps * x / pow((1 + x * x), 3.0 / 2.0);
    return tmp;

}
