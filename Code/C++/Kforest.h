#ifndef Kforest_H
#define Kforest_H

#include <iostream>
#include <cmath>
#include "forest.h"
#include "meteorology.h"

 class Kforest
 :
   public forest
 {
    const string name_;
    const string ofilename;
    double sigma_u;
    double sigma_v;
    double sigma_w;
    double L;
    double propConst;
    int timeStep_;
    int vertTimeStep_;
    double tau_u;
    double tau_v;
    double tau_w;
   public:
    Kforest(const dictionary&);
    void K ( double*, const double& );
    void GradK( double*, const double& );
    void updateParams ( ifstream* );
    void writeProfile ();
    inline void calcParams ( int, double* )
    {
	// not implemented; see Uforest class
    }

    inline string name()
    {
	return name_;
    }

    inline int horTimeStep()
    {
        return timeStep_;
    }

    inline int vertTimeStep()
    {
         return vertTimeStep_;
    }
 };

#endif
