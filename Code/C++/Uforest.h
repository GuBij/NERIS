#ifndef Uforest_H
#define Uforest_H

#include <iostream>
#include <cmath>
#include "forest.h"

 class Uforest
 :
   public forest
 {
    const string name_;
    const double pi=acos(-1.0);
    const double stackHeight;
    const string ofilename;
    const string Kfilename;
    const bool useLangevin;
    double dummy[4] ;
    double uStar;
    double mu;
    double zc;
    double zi;
    double windDirection;
    double elevation;
    double Uh;
    double GradUh;
    double intFGF( double );
    double heff_;
    double L;

   public:
    Uforest(const dictionary&);
    Uforest(const dictionary&, bool);
    void calcParams ( int, double* );
    void U ( double*, const double& );
    void updateParams ( ifstream* );
    void updateParams(double, double, double, double, double, double);
    void writeProfile ();
    inline string name()
    {
	return name_;
    }

    inline double wdir()
    {
        return windDirection;
    }

    inline double heff()
    {
        return heff_;
    }
 };

#endif
