#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

double ddot(double* v, double* w);
double* saxpy(double a, double* v, double* w);
double* matvec(double* v, int n);

double* cgsolve(int p, int n, int rank, double *b, int niters, int maxiterations )
{
    
double relres=1, alpha, beta,normb,rtrold,rtr;
double * Ad, * d, *x, *r;
double dot_product, dAd;
int i;
for(i=rank*n/p+1;i<=(rank+1)*n/p;i++)
{
    x[i-rank*n/p-1]=0;
    r[i-rank*n/p-1]=b[i-rank*n/p-1];
    d[i-rank*n/p-1]=r[i-rank*n/p-1];

}
dot_product = ddot(r,r);
MPI_Reduce(&dot_product, &rtr,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    if(rank==0)
{
    normb=sqrt(rtr);
}
while(relres > 1e-6 && niters < maxiterations)
{
    niters = niters+1;
    Ad = matvec(d,n);
    dot_product = ddot(d,Ad);
    MPI_Reduce(&dot_product,&dAd,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if(rank==0)
    {
        alpha = rtr / dAd;
        MPI_Bcast(&alpha,1, MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    x =saxpy(alpha,x,d);
    r =saxpy(-alpha,r,Ad);
    if(rank==0)
    {
        rtrold = rtr;
    }
    dot_product=ddot(r,r);
    MPI_Reduce(&dot_product,&rtr,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if(rank==0)
    {
        beta = rtr / rtrold;
        MPI_Bcast(&beta,1, MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    d = saxpy(beta,r,d);
    if(rank==0)
    {
        relres = sqrt(rtr) / normb;
        MPI_Bcast(&relres,1, MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
}

return x;
    
}

double ddot(double* v, double* w)
{
    if(sizeof(v)!=sizeof(w))
        printf("the size of two vectors does not match each other!\n");
    double sum=0;
    int i;
    for(i=0;i<=sizeof(v);i++)
        sum+=v[i]*w[i];
        return sum;
}

double* saxpy(double a, double* v, double* w)
{
    double* answer;
    int i;
    for(i=0;i<=sizeof(v);i++)
        answer[i]=v[i]+a*w[i];
    return answer;
}
double* matvec(double* v, int n)
{
    return v;
}

