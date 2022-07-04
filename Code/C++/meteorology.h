#ifndef meteorology_H
#define meteorology_H

#include <iostream>
#include "dictionary.h"

 class meteorology
 {
     string type;
   public:
     meteorology ();
     static meteorology* New(bool);
     static meteorology* New(string,bool);
     virtual string name() =0;
     virtual string meteofile() =0;
     virtual unsigned noRecords_() =0;
     virtual void U(double* u, const double& z)
     {
        u[0]=0; u[1]=0; u[2]=0;
     }
     virtual void K(double* k, const double& z)
     {
        k[0]=0; k[1]=0; k[2]=0;
     }
     virtual void GradK(double* k, const double& z)
     {
	k[0]=0; k[1]=0; k[2]=0;
     }

     virtual void updateTimeStep(double newTimeStep, double newTimeStepVert)
     {
	cerr << "Error: function 'updateParams ( ifstream* )' not implemented in class " << this->name() << endl;
     }

     virtual inline void corrCoeff(double* coeff)
     {
	coeff[0] = 0; coeff[1] = 0; coeff[2] = 0;
     }

     virtual inline double GradSig2(const double& z)
     {
        return 0.0;
     }

     virtual void updateParams(ifstream* ifile)
     {
	cerr << "Error: function 'updateParams ( ifstream* )' not implemented in class " << this->name() << endl;
     }

     virtual void calcParams(int timeIndex,double* data)
     {
	cerr << "Error: function 'updateParams ( ifstream* )' not implemented in class " << this->name() << endl;
     }

     virtual void updateParams (int timeIndex)
     {
	cerr << "Error: function 'updateParams ( int )' not implemented in class " << this->name() << endl;
     }

     virtual double LagrTimeScale (double z)
     {
	cerr << "Error: function 'LagrtimeScale ( double )' not implemented in class " << this->name() << endl;
     }

     virtual double wdir()
     {
	cerr << "Error: function 'wdir ( )' not implemented in class " << this->name() << endl;
     }

     virtual double heff()
     {
        cerr << "Error: function 'heff ( )' not implemented in class " << this->name() << endl;
     }

     virtual int horTimeStep()
     {
        cerr << "Error: function 'horTimeStep ( )' not implemented in class " << this->name() << endl;
     }

     virtual int vertTimeStep()
     {
        cerr << "Error: function 'vertTimeStep ( )' not implemented in class " << this->name() << endl;
     }
 };

#endif
