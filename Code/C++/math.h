#ifndef MATH_H
#define MATH_H

 #include <functional>
 using namespace std;

 namespace math
 {
   double residual
   (
	double* x,
        int nx, 
	double* F
   );

   int findZero
   (
	double* x,
        const int nx,
	function<void(double*,double*)> Jacobian,
	function<void(double*,double*)> Func 
   );

   int solveLinProblem
   (
        double* A,
        double* b,
        double* x,
        const int nx
   );

   int solveLinSystem
   (
	double* A, 
	double* b, 
	double* x,
        const int n
   );
 }

#endif
