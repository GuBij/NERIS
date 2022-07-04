#include <iostream>
#include <cmath>
#include "math.h"

namespace math
{
   const double f_SMALL_FLOAT = pow(2,-64);
   double max_noit = 10000;

   double residual(double* x, int nx, double* F)
   {
     double norm_func = 0, norm_x=0;
     for (int i = 0; i < nx; i++)
     {
       if ( abs(F[i]) > norm_func)
         norm_func = abs(F[i]);

       if ( abs(x[i]) > norm_x)
         norm_x = abs(x[i]);
     }
	return norm_func/norm_x;
   };

   int findZero(double* x, const int nx, function<void(double*,double*)> Jacobian, function<void(double*,double*)> Func)
   {
    int status=0;
    double dx[nx];
    double J[nx*nx];
    double F[nx];
    int noit = 0;

    Jacobian(x,J);
    Func(x,F);
    double fres = residual(x,nx,F);
    while (status == 0 && fres > f_SMALL_FLOAT && noit < max_noit)
    {
      status=solveLinProblem(J,F,dx,nx);
      for (int i=0; i < nx; i++)
      {
 	x[i] += dx[i];
      }
      Jacobian(x,J);
      Func(x,F);
      fres = residual(x,nx,F);
      noit++;
    }
    cout << "fres: " << fres << ", noit: " << noit << ", func: " << F[0] << ", jac: " << J[0] << "\n" << endl;
    return status;
  }

 int solveLinProblem(double* A, double* b, double* x, const int nx)
 {
    if (nx > 1)
	return solveLinSystem(A,b,x,nx);

    if (abs(A[0]) < f_SMALL_FLOAT)
	return 1;

    x[0] = b[0]/A[0];

    return 0;
 }

 int solveLinSystem(double* A, double* b, double* x, const int n) 
 {
    int rowIndex[n]; 
    const int noel = n*(n-1)/2;
    double rowEchelonMat[noel];
    for (int i = 0; i < n; i++)
    {
	rowIndex[i]=i;
    }

    for (int j = 0; j < (n-1); j++)
    {
	int rowIndexMax=rowIndex[j], swapIndex=j;
	for (int i = j+1; i < n; i++)
	{
	   if ( abs(A[j+rowIndex[i]*n]) > abs(A[j+rowIndexMax*n]) )
	   {
		rowIndexMax=rowIndex[i]; swapIndex=i;
	   }
	}
	rowIndex[swapIndex]=rowIndex[j]; rowIndex[j]=rowIndexMax;

	if ( abs(A[j+rowIndex[j]*n]) < f_SMALL_FLOAT )
	   return 1;

	for (int k = 1; k < n-j; k++)
	{
	  rowEchelonMat[k-1+noel-(n-j)*(n-j-1)/2]=A[k+j+rowIndex[j]*n]/A[j+rowIndex[j]*n];
	}

	for (int i = j+1; i < n; i++)
	{
	  double f = A[j+rowIndex[i]*n]/A[j+rowIndex[j]*n];

	  for (int k = j+1; k < n; k++)
	  {
	     A[k+rowIndex[i]*n] -= A[k+rowIndex[j]*n]*f;
	  }
	  b[rowIndex[i]] -= b[rowIndex[j]]*f;
	}
	b[rowIndex[j]] /= A[j+rowIndex[j]*n];
    }

    if ( abs(A[n-1+rowIndex[n-1]*n]) < f_SMALL_FLOAT )
	return 1;

    x[n-1]=b[rowIndex[n-1]]/A[n-1+rowIndex[n-1]*n];
    for (int i = 1; i < n; i++)
    {
	x[n-1-i]=b[rowIndex[n-1-i]];
	for (int j = 1; j < (i+1); j++)
	  x[n-1-i] -= rowEchelonMat[noel-i*(i+1)/2+j-1]*x[n-1-i+j];
    }

    return 0;
 };
}
