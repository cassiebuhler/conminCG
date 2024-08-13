/*********************************************************************/
/***    Copyright (c) Robert J. Vanderbei, 1994                    ***/
/***    All Rights Reserved                                        ***/
/*********************************************************************/

/* This is a function to return the elasped CPU time since this process
 * started running.  
 */

#ifdef MacOSX
#include <sys/time.h>
#else
#include <time.h>
#endif
#define CLOCKS_PER_SEC  1000000	/* [XSI] */

#ifdef QuadPrec
#include "Quad/Quad.h"
#define double Quad
#else
#define high(x) (x)
#endif

#ifdef __cplusplus
extern "C" {
#endif

double cputimer(
)
{
  double elaptime;

  elaptime = ((double)clock()) / CLOCKS_PER_SEC;

  return elaptime;
}

#ifdef __cplusplus
}
#endif
