#ifndef KTaylorOF_H
#define	KTaylorOF_H

 #include <iostream>
 #include <cmath>
 #include "openField.h"

 class KTaylorOF
 :
 public openField
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
   KTaylorOF(const dictionary&);
   void K( double*, const double& );
   void updateParams( ifstream* );
   double GradSig2( const double& );

   inline void updateTimeStep(double newTimeStep, double newTimeStepVert)
   {
	timeStep = newTimeStep; 
	timeStepVert = newTimeStepVert;
   }

   inline void corrCoeff(double* coeff)
   {
      coeff[0] = autoCorr[0]; coeff[1] = autoCorr[1]; coeff[2] = autoCorr[2];

      if ( abs(autoCorr[2]) > 1000 )
      {
	cout << "\ntau_w: " << tau_w << ", timeStepVert: " << timeStepVert << endl;
      }
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
