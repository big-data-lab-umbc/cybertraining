#include "memory.h"

/* 09/12/02-10/10/02, updated 02/07/08, 02/23/16 by Matthias K. Gobbert */

int *allocate_int_vector (int n)
{
  int *x;

  x = (int*) calloc (n, sizeof(int));

  if (x == NULL)
  {
    fprintf (stderr, "Problem allocating memory for vector\n");
#ifdef PARALLEL
    MPI_Abort (MPI_COMM_WORLD, 1);
#else
    exit (1);
#endif
  }

  return x;
}


long int *allocate_long_int_vector (long int n)
{
  long int *x;

  x = (long int*) calloc (n, sizeof(long int));

  if (x == NULL)
  {
    fprintf (stderr, "Problem allocating memory for vector\n");
#ifdef PARALLEL
    MPI_Abort (MPI_COMM_WORLD, 1);
#else
    exit (1);
#endif
  }

  return x;
}

double *allocate_double_vector (long n)
{
  double *x;

  x = (double*) calloc (n, sizeof(double));

  if (x == NULL)
  {
    fprintf (stderr, "Problem allocating memory for vector\n");
#ifdef PARALLEL
    MPI_Abort (MPI_COMM_WORLD, 1);
#else
    exit (1);
#endif
  }

  return x;
}

void free_vector (void *x)
{
  free (x);
}
