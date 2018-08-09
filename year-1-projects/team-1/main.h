#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "memory.h"
#include "cg.h"

typedef struct
{
  double tau, D, xmax, xmin, L;   
} S ;

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifndef M_PI
#define M_PI          (3.1415926535897932384626433832795029)
#endif

double frhs (double x, double y, double t, double tau, double D, double L);
double utrue (double x, double y, double t, double tau, double D, double L);

#endif

