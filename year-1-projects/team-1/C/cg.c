/* This function is designed for the Poisson problem discretized by
   a mesh with N points in each spatial dimension and a total of
   n mesh points in the mesh, that is, n = N^2 in two and n = N^3
   in three dimensions, for instance. The problem is assumed to be
   split in the last spatial dimension and l_N = N / np and l_n = n / np.
   The inputs to this function include l_n, l_N, N, and np,
   which are assumed to be supplied consistently.
   The vectors l_x, l_r, l_p, l_q, gl, gr should be defined and allocated
   in the calling routine. The vectors l_x, l_r, l_p, l_q should be of
   length l_n and gl and gr should be of length l_n / l_N.

   This function solves Ax = b for x using the (unpreconditioned)
   conjugatre gradient (CG) method. The matrix-vector product y = A * x
   for vectors x and y stored in l_x and l_y on each MPI process
   must be supplied in a function Ax() in a file Ax.c with header file Ax.h
   and have the prototype

   void Ax(double *l_y, double *l_x, long l_n, long l_N, long N,
     int id, int idleft, int idright, int np, MPI_Comm comm,
     double *gl, double *gr);

   The user must also supply a parallel dot product with prototype

   double parallel_dot(double *l_x, double *l_y, long l_n, MPI_Comm comm);

   We solve for l_x. The initial guess is supplied to the function as the 
   initial content of l_x. flag, relres, and iter are output variables that
   describe the result of the call. l_r is the right-hand side l_b
   when the function is called and contains output data of l_r
   when the call is complete .

   This implementation of the CG algorithm uses only 4 vectors:
   l_x = initial guess on input, solution on output,
   l_r = right-hand side l_b on input, residual on output,
   l_p = auxiliary variable = search direction,
   l_q = auxiliary variable = local portion of q = A * p
   These vectors are of length l_n.
   Two auxiliary vectors gl and gr of length l_n / l_N must be supplied. */

#include "cg.h"

void cg(double *l_x, int *flag, double *relres, int *iter, /* output */
        double *l_r, double tol, int maxit, /* input */
        double *l_p, double *l_q,
        long l_n, long l_N, long N,
        int id, int idleft, int idright, int np, MPI_Comm comm,
        double *gl, double *gr, double dt, double h, double D) {
  
  int it;
  long l_i;
  double n2b, tolb, normr, alpha, pq, beta, rho, rho1;

  n2b = sqrt(parallel_dot(l_r,l_r,l_n,comm)); /* Norm of the rhs vector b */
  if(n2b <= 1.0e-14) {           /* if the rhs vector is all zeros...     */
    for (l_i=0; l_i<l_n; l_i++)  /* then the solution is all zeros        */
      l_x[l_i] = 0.0;
    *flag = 0;                   /* a valid solution has been obtained    */
    *relres = 0;                 /* the relative residual is actually 0/0 */
    *iter = 0;                   /* no iterations need to be performed    */
    return;
  }

  tolb = tol * n2b; /* relative tolerance */
  Ax(l_q, l_x, l_n,l_N,N,id,idleft,idright,np,comm,gl,gr,dt,h,D); /* q = A * x   */
  for (l_i=0; l_i<l_n; l_i++)                        /* r = r - A*x       */
    l_r[l_i] = l_r[l_i] - l_q[l_i];
  rho = parallel_dot(l_r, l_r, l_n, comm);           /* rho = r' * r      */
  normr = sqrt(rho);                                 /* normr = norm(r)   */
  
  /* initial guess is a good enough: */
  if(normr <= tolb) {
    *flag = 0;
    *relres = normr / n2b;
    *iter = 0;
    return;
  }

  it = 0;
  while ((normr > tolb) && (it < maxit)) {
    it++;

    if(it == 1) {
#ifdef BLAS
      cblas_dcopy(l_n,l_r,1,l_p,1);                  /* p = r             */
#else
      // #pragma omp parallel for
      for (l_i=0; l_i<l_n; l_i++)                    /* p = r             */
        l_p[l_i] = l_r[l_i];
#endif
    } else {
      beta = rho / rho1;
#ifdef BLAS
      cblas_daxpby(l_n,1.0,l_r,1,beta,l_p,1);        /* p = r + beta * p  */
#else
      // #pragma omp parallel for
      for (l_i=0; l_i<l_n; l_i++)                    /* p = r + beta * p  */
        l_p[l_i] = l_r[l_i] + beta * l_p[l_i];
#endif
    }

    Ax(l_q, l_p, l_n,l_N,N,id,idleft,idright,np,comm,gl,gr,dt,h,D); /* q = A * p */

    pq = parallel_dot(l_p, l_q, l_n, comm);          /* pq = p' * q;      */
    alpha = rho / pq;                                /* alpha = rho/p'*q  */

#ifdef BLAS
    cblas_daxpy(l_n,(-alpha),l_q,1,l_r,1);           /* r = r - alpha * q */
#else
    // #pragma omp parallel for
    for (l_i=0; l_i<l_n; l_i++)                      /* r = r - alpha * q */
      l_r[l_i] = l_r[l_i] - alpha * l_q[l_i];
#endif
    rho1 = rho;                                      /* rho1 = rho;       */
    rho = parallel_dot(l_r, l_r, l_n, comm);         /* rho = r' * r */
    normr = sqrt(rho);                               /* normr = norm(r)   */

#ifdef BLAS
    cblas_daxpy(l_n,alpha,l_p,1,l_x,1);              /* x = x + alpha * p */
#else
    // #pragma omp parallel for
    for (l_i=0; l_i<l_n; l_i++)                      /* x = x + alpha * p */
      l_x[l_i] = l_x[l_i] + alpha * l_p[l_i];
#endif

    /* if (id == 0)
      printf("%d: normr <=? tolb: %12.7e %12.7e\n", it, normr, tolb); */
  }

  /* when the loop completes it will break to here - set output variables */  
  if (it < maxit)
    *flag = 0;
  else
    *flag = 1; /* maxit iterations reached */
  *relres = normr / n2b;
  *iter = it;

}

