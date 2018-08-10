#ifndef AX_H
#define AX_H

#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "memory.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

void Ax(double *y, double *x, long l_n, long l_N, long N,
     int id, int idleft, int idright, int np, MPI_Comm comm,
     double *gl, double *gr, double dt, double h, double D);

#endif

