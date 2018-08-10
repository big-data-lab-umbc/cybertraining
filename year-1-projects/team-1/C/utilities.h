#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

double serial_dot (double *x, double *y, long n);
double parallel_dot(double *l_x, double *l_y, long l_n, MPI_Comm comm);

#endif

