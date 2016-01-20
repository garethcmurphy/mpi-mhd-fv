/* $Id: emf_exchange.cpp,v 1.5 2006-10-30 15:17:55 gmurphy Exp $  */

#include "mpi.h"
#include "out.h"
#include "problem.h"
#include "physics.h"
#define MWULD MPI_COMM_WORLD
int
emf_exchange (TNT::Array2D < double >faceEx,
	      TNT::Array2D < double >faceEy,
	      int myNorth, int mySouth, int myEast, int myWest, int myid)
{


	int nx=faceEx.dim1();
	int ny=faceEy.dim2();
  MPI_Status stat;

  int fnx = faceEx.dim1 ();
  int fny = faceEy.dim2 ();

  int xfnx = faceEx.dim1 ();
  int xfny = faceEx.dim2 ();

  int yfnx = faceEy.dim1 ();
  int yfny = faceEy.dim2 ();

  int ne = NE;

  unk test;


  double gammam1 = PhysConsts::gamma - 1;
  double gammam1i = 1.0 / gammam1;


  MPI_Datatype FaceCol;
  MPI_Datatype FaceRowt;

  MPI_Datatype xFaceCol;
  MPI_Datatype xFaceRow;
  MPI_Datatype yFaceCol;
  MPI_Datatype yFaceRow;








  MPI_Type_vector (fnx,		/* # column elements */
		   1,		/* 1 column only */
		   fny,		/* skip ny elements */
		   MPI_DOUBLE,	/* elements are double */
		   &FaceCol);	/* MPI derived datatype */

  MPI_Type_vector (1,		/* # row elements */
		   fny,		/* 1 row only */
		   0,		/* skip ny elements */
		   //    myrow,   /* elements are double */
		   MPI_DOUBLE,	/* elements are double */
		   &FaceRowt);	/* MPI derived datatype */



  MPI_Type_vector (xfnx,	/* # column elements */
		   1,		/* 1 column only */
		   xfny,	/* skip ny elements */
		   MPI_DOUBLE,	/* elements are double */
		   &xFaceCol);	/* MPI derived datatype */

  MPI_Type_vector (1,		/* # row elements */
		   xfny,	/* 1 row only */
		   0,		/* skip ny elements */
		   //    myrow,   /* elements are double */
		   MPI_DOUBLE,	/* elements are double */
		   &xFaceRow);	/* MPI derived datatype */



  MPI_Type_vector (yfnx,	/* # column elements */
		   1,		/* 1 column only */
		   yfny,	/* skip ny elements */
		   MPI_DOUBLE,	/* elements are double */
		   &yFaceCol);	/* MPI derived datatype */

  MPI_Type_vector (1,		/* # row elements */
		   yfny,	/* 1 row only */
		   0,		/* skip ny elements */
		   //    myrow,   /* elements are double */
		   MPI_DOUBLE,	/* elements are double */
		   &yFaceRow);	/* MPI derived datatype */


  MPI_Type_commit (&FaceCol);
  MPI_Type_commit (&FaceRowt);
  MPI_Type_commit (&yFaceCol);
  MPI_Type_commit (&yFaceRow);
  MPI_Type_commit (&xFaceCol);
  MPI_Type_commit (&xFaceRow);



#define DEBUG
#ifdef DEBUG1
  std::cout << "Proc " << myid << " okay to here " << std::endl;
#endif


  //========== FACE CENTRED ELECTRIC FIELD EXCHANGE ===============

  // Send Columns

#define OLDWAY
#ifdef OLDWAY

  MPI_Send (&faceEx[0][xfny - 3], 1, xFaceCol, myNorth, myid, MWULD);
  MPI_Recv (&faceEx[0][1], 1, xFaceCol, mySouth, mySouth, MWULD, &stat);
  MPI_Barrier (MWULD);
  MPI_Send (&faceEx[0][2], 1, xFaceCol, mySouth, myid, MWULD);
  MPI_Recv (&faceEx[0][xfny - 2], 1, xFaceCol, myNorth, myNorth, MWULD, &stat);
  MPI_Barrier (MWULD);
  MPI_Send (&faceEy[2][0], 1, yFaceRow, myWest, myid, MWULD);
  MPI_Recv (&faceEy[yfnx - 2][0], 1, yFaceRow, myEast, myEast, MWULD, &stat);
  MPI_Barrier (MWULD);
  MPI_Send (&faceEy[yfnx - 3][0], 1, yFaceRow, myEast, myid, MWULD);
  MPI_Recv (&faceEy[1][0], 1, yFaceRow, myWest, myWest, MWULD, &stat);
  MPI_Barrier (MWULD);


#else


#ifdef PERIODIC
  for (int i = 2; i < nx-1 ; ++i)
  {
     faceEx[i][ny-2]=faceEx[i][2] ;
     faceEx[i][1]= (faceEx[i][ny-3]) ;

  }
  for (int j = 2; j < ny - 1; ++j)
  {
     faceEy[nx-2][j]= faceEy[2][j];
     faceEy[1][j]= (faceEy[nx-3][j]) ;
  }
#else

  for (int i = 2; i < nx-1 ; ++i)
  {
     faceEx[i][ny-2]=faceEx[i][ny-3] ;
     faceEx[i][1]= (faceEx[i][2]) ;
  }
  for (int j = 2; j < ny - 1; ++j)
  {
     faceEy[nx-2][j]= faceEy[nx-3][j];
     faceEy[1][j]= (-faceEy[2][j]) ;
  }

#endif

#endif


#ifdef CYLINDRICAL


if (myNorth == MPI_PROC_NULL)
{
  for (int i = 0; i < nx ; ++i)
  {
     faceEx[i][ny-2]=faceEx[i][ny-3] ;
  }
}

if (mySouth == MPI_PROC_NULL)
{
  for (int i = 0; i < nx ; ++i)
  {
     faceEx[i][1]= (faceEx[i][2]) ;
  }
}

if (myEast == MPI_PROC_NULL)
{
  for (int j = 0; j < ny ; ++j)
  {
     faceEy[nx-2][j]= faceEy[nx-3][j];
  }
}


if (myWest == MPI_PROC_NULL)
{
  for (int j = 0; j < ny ; ++j)
  {
     faceEy[1][j]= (-faceEy[2][j]) ;
  }
}
#endif

#ifdef MAES


  if (mySouth == MPI_PROC_NULL)
    {
      if (myWest == MPI_PROC_NULL)
   {
		int diskheight=8;
		int boxwidth=8;
     for (int j = 0; j < diskheight; ++j)
       {
         for (int i = 0; i < boxwidth; ++i)
      {
        faceEx[i][diskheight-1] = faceEx[i][diskheight];
        faceEx[i][j] = faceEx[i][diskheight];
        faceEy[i][j] = faceEy[i][diskheight];
      }
       }
     for (int j = 0; j < diskheight; ++j)
       {
        faceEy[boxwidth-1][j] = faceEy[boxwidth][j];
		 }
   }
    }


#endif



  return 0;
}
