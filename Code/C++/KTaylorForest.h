#ifndef KTaylorForest_H
#define	KTaylorForest_H

 #include <iostream>
 #include <cmath>
 #include "forest.h"

 class KTaylorForest
 :
 public forest
 {
   const string name_;
   const string meteoreologyfile_;
   int timeStep;
   int timeStepVert;
   double tau_u;
   double tau_v;
   double tau_w;
   double autoCorr[3];
   double sd_u;
   double sd_v;
   double sd_w;
   double hABL;
   double propConst;
   double L;

  public:
   KTaylorForest( const dictionary& );
   void K( double*, const double& );
   double GradSig2( const double& );
   void updateParams( ifstream* );

   inline void updateTimeStep(double newTimeStep, double newTimeStepVert)
   {
	timeStep = newTimeStep; 
	timeStepVert = newTimeStepVert;
   }

   inline void corrCoeff(double* coeff)
   {
      coeff[0] = autoCorr[0]; coeff[1] = autoCorr[1]; coeff[2] = autoCorr[2];
   }

   inline string name()
   {
       return name_;
   }

   inline int horTimeStep()
   {
        return timeStep;
   }

   inline int vertTimeStep()
   {
        return timeStepVert;
   }

 };

#endif
