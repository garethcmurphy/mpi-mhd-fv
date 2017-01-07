#include "TwoDim.h"

double **
TwoDim(int nrow, int ncol) {
    double **m;
    int i;

    m = (double **) malloc((size_t) (nrow * sizeof(double *)));
    m[0] = (double *) malloc((size_t) ((nrow * ncol) * sizeof(double)));

    for (i = 1; i < nrow; i++)
        m[i] = m[i - 1] + ncol;

    return m;
}
