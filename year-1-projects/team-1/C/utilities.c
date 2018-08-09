#include "utilities.h"

/* 09/12/02-10/10/02, updated 02/21/08 and 02/25/09 by Matthias K. Gobbert */

double serial_dot (double *x, double *y, long n)
{
  double dp;
  long i;

  dp = 0.0;
  for (i=0;i<=n;i++)
    dp = dp + x[i]*y[i];    

  return dp;
}

double parallel_dot(double *l_x, double *l_y, long l_n, MPI_Comm comm) {

  double l_dot=0.0;
  double dot=0.0;

  l_dot=serial_dot(l_x,l_y,l_n);

  MPI_Allreduce(&l_dot, &dot, 1, MPI_DOUBLE, MPI_SUM, comm);

  return dot;
}
