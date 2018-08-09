#ifndef CG_H
#define CG_H

#include <math.h>
#include "Ax.h"
#include "utilities.h"

void cg(double *l_x, int *flag, double *relres, int *iter, /* output */
        double *l_r, double tol, int maxit, /* input */
        double *l_p, double *l_q,
        long l_n, long l_N, long N,
        int id, int idleft, int idright, int np, MPI_Comm comm,
        double *gl, double *gr, double dt, double h, double D);

#endif

