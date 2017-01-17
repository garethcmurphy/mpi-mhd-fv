

#include "riemann.h"
#include "gtest/gtest.h"


TEST(sgn_test, int_arr_sort
) {
    double num = -7;
    double result = sgn(num);
    EXPECT_EQ(result,
              -1);
}


TEST(riemann_test, int_arr_sort
) {
    double num = -7;
    double lhs[] = {1, 0, 0, 0, 1, 1, 1, 1};
    double rhs[] = {1, 0, 0, 0, 1, 1, 1, 1};
    double flux[7];
    double res_state[7];
    double exp_res_state[] = {0, 0, 0, 0, 0, 0, 0, 0};
    double *max_speed;
    int idir;
    int time_step;
    double result = riemann(lhs, rhs, flux, res_state, time_step, max_speed, idir);
    EXPECT_FLOAT_EQ(exp_res_state[0], res_state[0]
    );
}
