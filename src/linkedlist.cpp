/* $Id: linkedlist.cpp,v 1.3 2006-10-30 15:17:55 gmurphy Exp $  */
#include <iostream>
#include <fstream>
#include <sstream>
#include "main.h"
#include "falle.h"
using namespace std;

int ne = 0;
int nx = 512;
int ny = 0;
int nz = 0;
double delta_x = 1.0;

double gammag = 1.666666666666667;
double gammam1 = gammag - 1.;
double gammam1i = 1. / gammam1;

double max (double a, double b);



struct unk_zone
{
  struct unk_zone *next;
  struct unk_zone *prev;
  struct unk_zone *up;
  double array[8];
  int loc;
};

typedef struct unk_zone unk;

int copy_zone (unk * mid_cell, double *state, int loc);
int print_list (char *filename, unk * mesh, int step);
int assign (unk * mesh, int size);



int
main ()
{
  unk *mesh;
  unk *soln;
  unk *ptr;
  unk *left_cell;
  unk *mid_cell;
  unk *right_cell;
  unk *store_cell;
  int ii = 0;
  int rc = 0;

  double left[8];
  double right[8];
  double lflux[8];
  double rflux[8];
  double maxspeed = 0;
  double delta = 0;
  double rho, vx, vy, vz, p, bx, by, bz;
  double kinetic, magnetic, energy;
  double time = 0, delta_x = 1 / 512., delta_t = 0, cfl = 0.8;
  ofstream outFile;
  int step = 1;
  int counter = 0;

  store_cell = (unk *) malloc (sizeof (unk));
  left_cell = (unk *) malloc (sizeof (unk));
  right_cell = (unk *) malloc (sizeof (unk));
  ptr = (unk *) malloc (sizeof (unk));

  outFile.open ("output/falle.log");
  outFile.close ();


  rho = 1;
  vx = 10;
  vy = 0;
  vz = 0;
  p = 20;
  bx = 1.414;
  by = 1.414;
  bz = 0.;
  kinetic = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
  magnetic = 0.5 * (bx * bx + by * by + bz * bz);

  left[0] = rho;
  left[1] = rho * vx;
  left[2] = rho * vy;
  left[3] = rho * vz;
  left[4] = kinetic + p * gammam1i + magnetic;
  left[5] = bx;
  left[6] = by;
  left[7] = bz;


  rho = 1;
  vx = -10;
  vy = 0;
  vz = 0;
  p = 1;
  bx = 1.414;
  by = 1.414;
  bz = 0.;
  kinetic = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
  magnetic = 0.5 * (bx * bx + by * by + bz * bz);


  right[0] = rho;
  right[1] = rho * vx;
  right[2] = rho * vy;
  right[3] = rho * vz;
  right[4] = kinetic + p * gammam1i + magnetic;
  right[5] = bx;
  right[6] = by;
  right[7] = bz;

  mesh = (unk *) malloc (sizeof (unk));
  mesh->next = NULL;
  mesh->prev = NULL;
  mesh->loc = 0;
  rc = copy_zone (mesh, left, 0);

  rc = assign (soln, nx);
  // Assign a 50 element array of unks
  for (ii = 1; ii < nx; ii++)
    {
      // Allocate a temporary variable
      mid_cell = (unk *) malloc (sizeof (unk));
      mid_cell->next = NULL;
      mid_cell->loc = ii;
      if (ii < nx / 2.0)
	rc = copy_zone (mid_cell, left, ii);
      else
	rc = copy_zone (mid_cell, right, ii);
      for (ptr = mesh; ptr->next != NULL; ptr = ptr->next)
	;
      ptr->next = mid_cell;
      mid_cell->prev = ptr;
    }

  step = 0;
  rc = print_list ("out_00", mesh, step);

  for (step = 1; step < 900; step++)
    {
// Advect the linked list
//               counter=0;
      for (ptr = mesh; ptr->next != NULL; ptr = ptr->next)
	{
//              counter++;

	  rc = maxspd (ptr->array, &maxspeed);
	}
      // rc = print_list ("out_00", mesh, 100*step);


      cfl = 0.8;
      delta = cfl / maxspeed;
      delta_x = 1 / 512.0;
      delta_t = delta_x * delta;
      time = time + delta_t;
      cout
	<< "Step: " << step
	<< "\tTime: " << time
	<< "\tTstep: " << delta_t << "\tMax speed: " << maxspeed
//      << "\tCounter: " << counter 
	<< endl;

      //    mid_cell = mesh;
      //   right_cell = mesh;
      //  left_cell = mesh;
      // stre_cell is reinitialised here to mesh
      //   store_cell = mesh;
      //  rc = copy_zone (store_cell, left, 0);

      int idir = 1;
      // Loop through the mesh
      for (ptr = mesh; ptr->next != NULL; ptr = ptr->next)
	{
	  if (ptr->prev != NULL)
	    {
	      left_cell = ptr->prev;
	      mid_cell = ptr;
	      right_cell = ptr->next;

	      // Riemann solve between mid_cell and left_cell
	      rc =
		falle_mhd (left_cell->array, mid_cell->array, lflux, step,
			   &maxspeed, idir);

	      // Riemann solve between right_cell and left_cell
	      rc =
		falle_mhd (mid_cell->array, right_cell->array, rflux, step,
			   &maxspeed, idir);

	      if (left_cell->prev != NULL)
		{
		  rc =
		    copy_zone (left_cell, store_cell->array, store_cell->loc);
//                      left_cell->array[0]=999;
		}

	      for (ii = 0; ii < 8; ii++)
		{
		  store_cell->array[ii] = mid_cell->array[ii]
		    - delta * (rflux[ii] - lflux[ii]);
		}
	      store_cell->loc = mid_cell->loc;

	    }


	}
      if (step % 10 == 0)
	{
	  rc = print_list ("out_00", mesh, step);
	}
    }

// Print out the linked list
  outFile.open ("out.dat");
  for (left_cell = mesh; left_cell->next != NULL; left_cell = left_cell->next)
    {
      for (ii = 0; ii < 8; ii++)
	{
	  outFile << " " << left_cell->array[ii];
	}
      outFile << endl;
      cout << " " << left_cell->array[0];
    }
  outFile.close ();

  cout << endl;
  return 0;
}

int
copy_zone (unk * mid_cell, double *state, int loc)
{
  int ii;
  mid_cell->loc = loc;
  for (ii = 0; ii < 8; ii++)
    mid_cell->array[ii] = state[ii];
}

int
print_list (char *filename, unk * mesh, int step)
{
  int ii = 0;
  unk *left_cell;
  ofstream outFile;
  double rho, rhoi, vx, vy, vz, p, et, ke;

  stringstream s;
  stringstream stream_filename;
  stringstream stream_temp_b;
  string str_file_tag;
  string str_output_filename;
  string str_input_filename;
  string outputdir = "output/";


  s.clear ();
  s.width (5);
  s.fill ('0');
  s << step;
  s >> str_file_tag;
  stream_filename.clear ();
  stream_filename << outputdir << filename << str_file_tag;
  stream_filename >> str_input_filename;





  outFile.open (str_input_filename.c_str ());
  for (left_cell = mesh; left_cell->next != NULL; left_cell = left_cell->next)
    {
      rho = left_cell->array[0];
      rhoi = 1 / rho;
      vx = left_cell->array[1] * rhoi;
      vy = left_cell->array[2] * rhoi;
      vz = left_cell->array[3] * rhoi;
      outFile << " " << left_cell->loc;
      outFile
	<< " " << rho
	<< " " << vx
	<< " " << vy
	<< " " << vz
	<< " " << left_cell->array[4]
	<< " " << left_cell->array[5]
	<< " " << left_cell->array[6] << " " << left_cell->array[7];
      outFile << endl;
//      cout << " " << left_cell->array[0];
    }
  outFile.close ();
}

int
assign (unk * mesh, int size)
{

  int ii = 0;
  unk *temp;
  unk *ptr;

  mesh = (unk *) malloc (sizeof (unk));
  mesh->next = NULL;
  mesh->prev = NULL;
  for (ii = 1; ii < size; ii++)
    {
      // Allocate a temporary variable
      temp = (unk *) malloc (sizeof (unk));
      temp->next = NULL;
      temp->loc = ii;
      for (ptr = mesh; ptr->next != NULL; ptr = ptr->next)
	;
      ptr->next = temp;
      temp->prev = ptr;
    }
}
