#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>


#ifdef HAVE_LIBMPICH
#include "mpi.h"
#endif

#define TWODIM

#include "tnt.h"
using namespace TNT;

using namespace std;
#define N_X 100
#define N_Y 10
#define N_E 4
#define DEBUG 1
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
//#undef DEBUG
#ifndef GLOBAL_VARIABLES
#define GLOBAL_VARIABLES 1
extern int nx;
extern int ny;
extern int ne;
extern double delta_x;
extern double gammag;
extern double gammam1;
extern double gammam1i;


#define _MASS .array[0]
#define _MOMX .array[1]
#define _MOMY .array[2]
#define _MOMZ .array[3]
#define _ENER .array[4]
#define _B_X .array[5]
#define _B_Y .array[6]
#define _B_Z .array[7]
#define _COOLING .temperature

typedef struct
{
  double array[8];
  double temperature;
  double cooling;
} zone;




#endif //  GLOBAL_VARIABLES
