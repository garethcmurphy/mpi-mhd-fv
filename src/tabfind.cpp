#include "tabfind.h"


void
tabfind(double temperature, double *rate, double *atomic_temp_tab) {


    double cooling_tab[50] = {
            -22.97, -22.87, -22.77, -22.64, -22.52, -22.42, -22.31, -22.22, -22.14,
            -22.07, -22.02, -21.97, -21.91, -21.86, -21.81, -21.75, -21.69, -21.63,
            -21.58, -21.53, -21.49, -21.45, -21.41, -21.38, -21.37, -21.36, -21.37,
            -21.4, -21.45, -21.52, -21.61, -21.69, -21.75, -21.78, -21.81, -21.85,
            -21.91, -21.96, -21.98, -21.98, -21.98, -21.97, -21.95, -21.95, -21.99,
            -22.1, -22.23, -22.34, -22.41, -22.44
    };
    int j = 0;
    double temperature1, temperature2, L1, L2, del_temperature, delL;
    for (j = 0; j < 50; j++) {
        cooling_tab[j] = cooling_tab[j] + 22.0;    // removing 10^-22 factor
    }

    locate(temperature, atomic_temp_tab, &j, 50);    //  temperature is between j and j+1
//      write(*,*) 'tab j', j
    if (j == 50) {
        *rate = cooling_tab[j];
    } else {
        temperature1 = atomic_temp_tab[j];
        temperature2 = atomic_temp_tab[j + 1];
        del_temperature = temperature2 - temperature1;
        L1 = cooling_tab[j];
        L2 = cooling_tab[j + 1];
        delL = L2 - L1;
        *rate = L1 + (temperature - temperature1) * delL / del_temperature;
    }

    *rate = pow(10.0, *rate);    // Lambda_bar [10-22 already removed]
    *rate = 0.01 * *rate;        //absorb 10-2 factor to obtain LN
}
