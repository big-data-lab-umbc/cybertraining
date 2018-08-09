#include "main.h"

#define SQ(X) ((X)*(X))

int main (int argc, char *argv[]) {

  int id, np, namelen, idleft, idright;
  int rem, lenm;                    
  char name[MPI_MAX_PROCESSOR_NAME];
  char filename[64];
  char message[100];
  int iter, maxit, flag;
  long N, l_N, n, l_n, i, j, l_j, l_k;
  long l_ia, l_ib, nt, Nt;
  long l_Ntmp, l_ntmp;
  double *l_r, *l_p, *l_q, *l_un;
  double *x, *y;
  double *gl, *gr;
  double h, tol, relres, t, kappa, tau, tout, dtout, tnew, D, xmax, xmin, L;
  int ntout, cum_iter;
  double enorminf, l_enorminf, err_ij, tmp;
  double min_u, l_min_u, max_u, l_max_u, mass_t, l_mass_t;
  double start, end, dt, t_fin;
  FILE *idfil;
  MPI_Comm comm;
  MPI_Status status;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &id);
  MPI_Get_processor_name (name, &namelen);

  if (id > 0)
    idleft = id - 1;
  else
    idleft = MPI_PROC_NULL;
  if (id < np-1)
    idright = id + 1;
  else
    idright = MPI_PROC_NULL;
  comm = MPI_COMM_WORLD;
  
  /* process command-line inputs: */
  if (argc != 6){
    if (id == 0){
      printf ("Usage: ./pollution N0 tol maxit t_fin dt\n");
      printf (" with int N0, double tol, int maxit, double t_fin, double dt\n");
    }
    MPI_Abort (MPI_COMM_WORLD, 1);
  }
  N     = (long)(atof(argv[1])); /* TMP: rename N as N0 */
  tol   =       (atof(argv[2]));
  maxit =  (int)(atof(argv[3]));
  t_fin =       (atof(argv[4]));
  dt    =       (atof(argv[5]));

  /* Added new parameters*/    
  kappa = 10.0;
  tau = 8.0; /*Also eventually from command line*/
  D = 10.0;
  xmax = 50.0;
  xmin = -xmax;
  L = xmax - xmin;

  Nt = ceil(t_fin/dt);
  h = L/(N+1.0);  /*step size*/ /* TMP N.B.: still for Dirichlet BC */
  n = SQ(N);  /* n = N^2 */ /* TMP: rename n as N */
  
  /* in the following distribution, the first rem=N%np processes get
   * N/np+1 values, and the remaining np-rem ones get N/np,
   * for a total of rem*(N/np+1) + (np-rem)*(N/np)
   * = rem*(N/np) + rem + np*(N/np) - rem*(N/np)
   * = (rem+np-rem)*(N/np) + rem = np*(N/np) + rem = (N-rem) + rem = N */
  rem = N%np; /* remainder when dividing N values into np MPI processes */
  if (id < rem) { /* if id is 0, 1, ..., rem-1 then one more value: */
    l_N = N/np + 1;                  /* local number of values */
  } else {
    l_N = N/np;                      /* local number of values */
  }
  l_n = N * l_N;    /* size of local matrix l_u */

  /* Notice carefully: start and ending index are coded "C-style",
   * i.e., a for-loop should read: for(l_i=l_ia; l_i<l_ib; l_i++)
   * with a strictly less than comparison. */
  rem = N%np; /* remainder when dividing N values into np MPI processes */
  if (id < rem) { /* if id is 0, 1, ..., rem-1 then one more value: */
    l_ia = id*l_N;                   /* local starting index */
    l_ib = l_ia + l_N;               /* local ending index */
  } else {
    l_ia = rem*(l_N+1)+(id-rem)*l_N; /* local starting index */
    l_ib = l_ia + l_N;               /* local ending index */
  }
  /* test output: */
  if (np <= 16) { /* for more than np=16 this becomes unreadable */
    sprintf(message, "P%03d: l_n=%12ld; l_N=%6ld from l_ia=%6ld to l_ib=%6ld", id, l_n, l_N, l_ia, l_ib);
    if (id == 0) {
      for (i = 0; i < np; i++) {
        if (i > 0)
          MPI_Recv(message,100,MPI_CHAR,i,i,MPI_COMM_WORLD,&status);
        printf("[%3d] %s\n", id, message);
      }
    } else {
      lenm = 1 + strlen(message);
      MPI_Send(message,lenm,MPI_CHAR,0,id,MPI_COMM_WORLD);
    }
  }

  if (id == 0) {
    printf("pollution:\n");
    printf("kappa =%14.6e D =%14.6e tau =%14.6e\n", kappa, D, tau);
    printf("xmax =%14.6e t_fin =%14.6e\n", xmax, t_fin);
    printf("N =%6ld n =%12ld h =%14.6e\n", N, n, h);
    printf("dt =%14.6e Nt =%6ld\n", dt, Nt);
    printf("tol =%14.6e maxit =%5d\n", tol, maxit);
    printf("M_PI =%14.6e SQ(h) =%14.6e\n", M_PI, SQ(h));
  }

  l_r   = allocate_double_vector(l_n); /* right-hand side in, residual out */
  l_un  = allocate_double_vector(l_n);                                        /* Get vector for u_n*/
  l_p   = allocate_double_vector(l_n);
  l_q   = allocate_double_vector(l_n);
  x = allocate_double_vector(N);
  y = allocate_double_vector(N);
  gl = allocate_double_vector(n/N);
  gr = allocate_double_vector(n/N);

  for (i = 0; i < N; i++){
    x[i] = -xmax + h * (double)(i+1); /* TMP N.B.: still interior points only */
  }
  for (j = 0; j < N; j++){
    y[j] = -xmax + h * (double)(j+1);
  }

  /* initializations of the 4 large arrays in the code: */
  /*solution vector at timestep n*/
  for (l_j = 0; l_j < l_N; l_j++){
    for (i = 0; i < N; i++){
      l_k = i + N*l_j;
      l_un[l_k] = 0.0; /* N.B.: initial condition 0.0 hard-wired here */
    }
  }
  /* search direction l_p */
  for (l_j = 0; l_j < l_N; l_j++){
    for (i = 0; i < N; i++){
      l_k = i + N*l_j;
      l_p[l_k] = 0.0;
    }
  }
  /* auxiliary vector l_q */
  for (l_j = 0; l_j < l_N; l_j++){
    for (i = 0; i < N; i++){
      l_k = i + N*l_j;
      l_q[l_k] = 0.0;
    }
  }
  
  for (i = 0; i < N; i++) {
    gl[i] = 0.0;   /*Setting to zero temporary */
    gr[i] = 0.0;
  }

  if (id == 0) {
    printf("     n  t_n  it cumit       min(u)      max(u)     mass(u)    enorminf\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();  /* start time */

  /* Time loop starts here */
  tout = 4.0; /* first output time */
  dtout = 4.0; /* output time step thereafter */
  ntout = 1;
  cum_iter = 0;
  for (nt = 0; nt<Nt; nt++) {
    t = nt * dt;
    tnew = (nt + 1) * dt;

    /* residual vector l_r = right-hand side l_b = l_un + frhs(tnew): */
    for (l_j = 0; l_j < l_N; l_j++){
      j = l_j + l_ia; /* l_j + l_N*id; */
      for (i = 0; i < N; i++){
        l_k = i + N*l_j;
        l_r[l_k] = l_un[l_k] + dt * frhs(x[i],y[j],tnew, tau,D,L);
      }
    }

    cg(l_un, &flag, &relres, &iter, l_r, tol, maxit,
       l_p, l_q, l_n, l_N, N, id, idleft, idright, np, comm, gl, gr, dt, h, D);
    cum_iter += iter;

    /* output every dtout times, up to 6 times total: */
    if ( (tnew>=tout-1.0e-14) && (ntout<=6) ) {

      /* Calculate min_u, max_u, enorminf: */
      l_min_u = l_un[0];
      l_max_u = l_un[0];
      l_mass_t = 0.0; /* TMP: mass calculations not implemented yet */
      l_enorminf=0.0;
      for (l_j = 0; l_j < l_N; l_j++){
        j = l_j + l_ia;
        for (i = 0; i < N; i++){
          l_k=i+l_j*N;
          if (l_un[l_k] < l_min_u)
            l_min_u = l_un[l_k];
          if (l_un[l_k] > l_max_u)
            l_max_u = l_un[l_k];
          if (h == L / (N+1)) { /* mesh consists of interior points only */
            l_mass_t += 4.0 * l_un[l_k]; /* this assumes u=0 on boundary */
          } else if (h == L / (N-1)) { /* mesh contains the boundary points */
            /* TMP: the following does not give same result as Matlab! Why? */
            if ( (i>0) && (i<N) && (j>0) && (j<N) ) { /* interior point */
              l_mass_t += 4.0 * l_un[l_k];
            } else {
                     if ( (i==0) && (j==0) ) { /* corner point */
                l_mass_t += 1.0 * l_un[l_k];
              } else if ( (i==N) && (j==0) ) { /* corner point */
                l_mass_t += 1.0 * l_un[l_k];
              } else if ( (i==0) && (j==N) ) { /* corner point */
                l_mass_t += 1.0 * l_un[l_k];
              } else if ( (i==N) && (j==N) ) { /* corner point */
                l_mass_t += 1.0 * l_un[l_k];
              } else {                         /* edge point */
                l_mass_t += 2.0 * l_un[l_k];
              }
            }
          } else {
            printf("Error: This point should never by reached!\n");
          }
          // U True Value
          l_q[l_k] = utrue(x[i],y[j],tnew, tau,D,L);
          err_ij=fabs(l_q[l_k]-l_un[l_k]);
          if (err_ij>l_enorminf)
            l_enorminf=err_ij;
        }
      }
      MPI_Reduce(&l_min_u   ,&min_u   , 1,MPI_DOUBLE, MPI_MIN, 0, comm);
      MPI_Reduce(&l_max_u   ,&max_u   , 1,MPI_DOUBLE, MPI_MAX, 0, comm);
      MPI_Reduce(&l_mass_t  ,&mass_t  , 1,MPI_DOUBLE, MPI_SUM, 0, comm);
      MPI_Reduce(&l_enorminf,&enorminf, 1,MPI_DOUBLE, MPI_MAX, 0, comm);
      mass_t *= SQ(h/2.0);

      /* output one line of table: */
      if (id == 0) {
        printf("%6ld%5.1f%4d%6d%13.4e%12.4e%12.4e%12.4e\n",
               nt, tnew, iter, cum_iter, min_u, max_u, mass_t, enorminf);
      }

      /* increment output time step counter and output time step value: */
      ntout = ntout + 1;
      tout = tout + dtout;
    }

} //end of time loop
    
  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();  /* end time */

  if (id == 0) {
//    printf("N = %6ld, DOF = %12ld, ||u - u1|| = %24.16e, #iter =%d, wall clock time = %10.2f\n", N, n, enorminf,iter,end-start);
//    printf("t_fin = %24.16e, dt = %24.16e\n", t_fin, dt);
//    printf("flag = %1d, iter = %d\n", flag, iter);
//    printf("relres             = %24.16e\n", relres);
//      printf("h                  = %24.16e\n", h);
//      printf("h^2                = %24.16e\n", h*h);
//      printf("enorminf           = %24.16e\n", enorminf);
//    printf("C = enorminf / h^2 = %24.16e\n", (enorminf/(h*h)));
    printf("wall clock time    = %10.2f seconds\n", end-start);
    fflush(stdout);
  }

  MPI_Finalize();
  free_vector(l_q);
  free_vector(l_r);
  free_vector(l_p);
  free_vector(l_un);
  free_vector(x);
  free_vector(y);
  free_vector(gl);
  free_vector(gr);

  return 0;
}

double frhs (double x, double y, double t, double tau, double D, double L) {

  double f;
  /* Matlab code for f:
    F = ((2*t/S.tau^2)*exp(-(t/S.tau)^2)) ...
          * ((cos((pi/S.L)*X)).^2.*(cos((pi/S.L)*Y)).^2) ...
      - (S.D * (1-exp(-(t/S.tau)^2)) * (-2*(pi/S.L)^2)) ...
          * ( cos((2*pi/S.L)*X).*(cos((pi/S.L)*Y)).^2 ...
            + (cos((pi/S.L)*X)).^2.*cos((2*pi/S.L)*Y) );
  */
  f = ((2.0*t/SQ(tau))*exp(-SQ(t/tau)))
        * (SQ(cos((M_PI/L)*x)) * SQ(cos((M_PI/L)*y)))
    - (D * (1.0-exp(-SQ(t/tau))) * (-2.0*SQ(M_PI/L)))
       * (cos((2.0*M_PI/L)*x) * SQ(cos((M_PI/L)*y))
         + SQ(cos((M_PI/L)*x)) * cos((2.0*M_PI/L)*y) );
  return f;
}

double utrue (double x, double y, double t, double tau, double D, double L) {

  double u;
  /* Matlab code for u:
    U = (1-exp(-(t/S.tau)^2)) * ((cos((pi/S.L)*X)).^2.*(cos((pi/S.L)*Y)).^2);
  */
  u = (1.0-exp(-SQ(t/tau))) * (SQ(cos((M_PI/L)*x)) * SQ(cos((M_PI/L)*y)));
  return u;
}

