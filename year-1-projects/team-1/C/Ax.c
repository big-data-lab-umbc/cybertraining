#include "Ax.h"
#define SQ(X) ((X)*(X))


void Ax(double *l_v, double *l_u, long l_n, long l_N, long N,
        int id, int idleft, int idright, int np, MPI_Comm comm,
        double *gl, double *gr, double dt, double h, double D) {

  long i, l_j;
  double tmp;
  MPI_Status status;
  double Ddthh = D * (dt / SQ(h));

  for (l_j = 1; l_j < l_N-1; l_j++) {
    for (i = 0; i < N; i++) {
      tmp = (1.0 + (4.0*Ddthh)) * l_u[i+N*l_j];
      if (l_j >   0  ) tmp -= Ddthh * l_u[ i   +N*(l_j-1)];
      if (  i >   0  ) tmp -= Ddthh * l_u[(i-1)+N* l_j   ];
      if (  i <   N-1) tmp -= Ddthh * l_u[(i+1)+N* l_j   ];
      if (l_j < l_N-1) tmp -= Ddthh * l_u[ i   +N*(l_j+1)];
      l_v[i+N*l_j] = tmp;
    }
  }

  if((id%2)==0) { /* even-numbered processes */
    MPI_Recv(gl,N,MPI_DOUBLE,idleft,0,comm,&status);
    MPI_Recv(gr,N,MPI_DOUBLE,idright,0,comm,&status);
    MPI_Send(&(l_u[0]),N,MPI_DOUBLE,idleft,0,comm);
    MPI_Send(&(l_u[N*(l_N-1)]),N,MPI_DOUBLE,idright,0,comm);
  } else { /* odd-numbered processes */
    MPI_Send(&(l_u[0]),N,MPI_DOUBLE,idleft,0,comm);
    MPI_Send(&(l_u[N*(l_N-1)]),N,MPI_DOUBLE,idright,0,comm);
    MPI_Recv(gl,N,MPI_DOUBLE,idleft,0,comm,&status);
    MPI_Recv(gr,N,MPI_DOUBLE,idright,0,comm,&status);
  }  

  l_j=0;
    for (i = 0; i < N; i++) {
      tmp = (1.0 + (4.0*Ddthh)) * l_u[i+N*l_j];
      tmp -= (Ddthh) * gl[i];
      if (  i >   0  ) tmp -= Ddthh * l_u[(i-1)+N* l_j   ];
      if (  i <   N-1) tmp -= Ddthh * l_u[(i+1)+N* l_j   ];
      if (l_j < l_N-1) tmp -= Ddthh * l_u[ i   +N*(l_j+1)];
      l_v[i+N*l_j] = tmp;
    }

  l_j=l_N-1; 
    for (i = 0; i < N; i++) {
      tmp = (1.0 + (4.0*Ddthh)) * l_u[i+N*l_j];
      if (l_j >   0  ) tmp -= Ddthh * l_u[ i   +N*(l_j-1)];
      if (  i >   0  ) tmp -= Ddthh * l_u[(i-1)+N* l_j   ];
      if (  i <   N-1) tmp -= Ddthh * l_u[(i+1)+N* l_j   ];
      tmp -= (Ddthh) * gr[i];
      l_v[i+N*l_j] = tmp;
    }

}

