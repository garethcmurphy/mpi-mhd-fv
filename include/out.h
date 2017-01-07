/* $Id: out.h,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include <mpi.h>
#include "tnt.h"
#include "PhysConsts.h"
#include <iostream>
#define SMALLX 2
#define SMALLY 2
#define BIGX 4
#define BIGY 4

#define NE 10

#ifndef UNK
#define UNK
typedef struct
{
  int temperature;
  double array[NE];
} unk;

#define _MASS .array[0]
#define _MOMX .array[1]
#define _MOMY .array[2]
#define _MOMZ .array[3]
#define _ENER .array[4]
#define _B_X .array[5]
#define _B_Y .array[6]
#define _B_Z .array[7]
#define _DIVB .array[8]
#define _COOL .array[9]



struct flux_struct
{
  double array[NE];
};
typedef struct flux_struct flux;


#endif
//typedef struct unk_zone unk;
int out (int numprocs,
	 int myid,
	 MPI_Comm new_comm,
	 int ndims,
	 int *dim_size,
	 TNT::Array2D < unk > smallarr, int *arr_size, int step);
int outhdf5 (int numprocs,
	 int myid,
	 MPI_Comm new_comm,
	 int ndims,
	 int *dim_size,
	 TNT::Array2D < unk > smallarr, int *arr_size, int step);
#define NGC 2
