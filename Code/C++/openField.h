#ifndef openField_H
#define	openField_H

 #include <fstream>
 #include <cmath>
 #include "dictionary.h"
 #include "meteorology.h"

 class openField
 :
 public dictionary, public meteorology
 {
   const string name_;
   const string meteofile_;
   const double g;
   const double zref;
   const double zref_T;
   const double zref_T0;
   const double RH;
   const double pressure;
   const double rho;
   const double R;
   const double cp;
   const double d;
   const double flowRate;
   const double Ts;
   const bool useLangevin;
   unsigned noRecords;
   void picard(double, double, double, double&, double&, double&);
   void plumeRise(double, double, double, double, double&);
   void LagrangianTimeScale(double, double, double, double, double, double, double, double*, double, double&);
   void timeStepSelection(const int, const double*, const double, int*);

   protected:
    const double kappa;
    const double z0;
    const double Sct;
    const double stackHeight;
    double stabKernel(double, string, bool);
 
   public:
    openField(const dictionary&);
    openField(const dictionary&, bool);
    void calcParams(int,double*);
    inline void updateParams(ifstream*)
    {
	// not implemented;
    }

    inline virtual string name()
    {
       return name_;
    }

    inline string meteofile()
    {
       return meteofile_;
    }

    inline unsigned noRecords_()
    {
	return noRecords+1;
    }
 };

#endif
