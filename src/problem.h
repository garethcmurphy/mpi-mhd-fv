/* $Id: problem.h,v 1.5 2006-11-16 13:48:07 gmurphy Exp $  */

//#define BLAST
//#define ORSZAGTANG
//#define JET
#define MAES



#ifdef BLAST
#undef CYLINDRICAL
//#define ROE
#define HLLD
#define PERIODIC
#define LIMITER vanleer
//#define LIMITER minmod
#endif


#ifdef ORSZAGTANG
//#define ROE
#define PERIODIC
#undef CYLINDRICAL
//#define LIMITER minmod
#define HLLD
#define LAXFRIEDRICHS
#define LIMITER vanleer
#endif

#ifdef JET
//#define RESISTIVE
//#define GRAVITY
#define HLLD
//#define HLLD
#define CYLINDRICAL
#define LIMITER vanleer
//#define COOLING
#endif


#ifdef MAES
#define RESISTIVE
#define GRAVITY
#define CYLINDRICAL
//#define ROE
#define HLLD
//#define LAXFRIEDRICHS
#define LIMITER minmod
#endif
