#ifndef forest_H
#define forest_H

 #include <iostream>
 #include <fstream>
 #include <cmath>
 #include <gsl/gsl_sf_hyperg.h>
 #include "dictionary.h"
 #include "meteorology.h"
 #include "math.h"

 class forest
 :
 public dictionary, public meteorology
 {
   const string name_;
   const string meteofile_;
   const double alpha;
   const double stackHeight;
   const double diam;
   const double flowRate;
   const double Ts;
   double sparsityCorr;
   double muM;
   double muH;
   double uStar;
   double thetaStar;
   double L;
   double zc;
   double zi;
   double LADTop;
   double Uh;
   double GradUh;
   double Z0;
   double d;
   unsigned noRecords;

   protected:
    const double plantCover;
    const double zStar;
    const double canopyDepth;
    const double a2;
    const double a3;
    const double a4;
    const double a1;
    const double kappa;
    const double htree;
    const double TAD;
    const double z0;
    const double z0g;
    const double zref;
    const double g;
    const double RH;
    const double pressure;
    const double rho;
    const double R;
    const double cp;
    const double zref_T;
    const double zref_T0;
    void timeStepSelection(const int, const double*, const double, int*);

    inline double d_()
    {
	return d;
    }

   public:
    forest(const dictionary& dict); 
    double stabKernel(double, string, string);
    double stabKernel(double, string);
    void stateEqn(double*, double*);
    void JacStateEqn(double*, double*);
    double findMu(double);
    double findMu(double, string, double);
    double picard(double, double, double, double*);
    int calibrate(double, double, double, function<void(double*,double*)>, function<void(double*,double*)>, function<void(double*,double*)>, function<void(double*,double*)>);
    void plumeRise(double, double, double, double, double, double&);
    void LagrangianTimeScale(double, double, double, double, double, double, double, double*, double, double&);
    void printParams();

    inline void LAD2(double* y, double* F, const double* rhs)
    {
	if (y[0]>a4)
	  y[0] = 0.99*a4;
	else if (y[0]<0)
	  y[0] = 0.001;

	F[0] = -(a1*pow(y[0],a2)*pow(a4-y[0],a3) - *rhs)*(a1*pow(y[0],a2)*pow(a4-y[0],a3) - *rhs);
    }

    inline void JacLAD2(double* y, double* F, const double* rhs)
    {
	 if (y[0]>a4)
	  y[0] = 0.99*a4;
        else if (y[0]<0)
          y[0] = 0.001;

         F[0] = 2*(a1*pow(y[0],a2)*pow(a4-y[0],a3) - *rhs)*a1*(a2*pow(y[0],a2-1)*pow(a4-y[0],a3)-a3*pow(y[0],a2)*pow(a4-y[0],a3-1));
    }

    inline double roughSublayerCorr(double z, double mu_)
    {
	return 2.0/3.0*log(1+1.5*(zStar-d)/(mu_*z))*exp(-mu_*z/(zStar-d));
    }

    inline double roughSublayerCorr(double z, double mu_, string type, double L)
    {
          return stabKernel((1+0.5*(zStar-d)/(mu_*z))*z/L,type)*2.0/3.0*log(1+1.5*(zStar-d)/(mu_*z))*exp(-mu_*z/(zStar-d));
    }

    inline double JacRoughSublayerCorr(double z, double mu_)
    {
	return -exp(-mu_*z/(zStar-d))*((zStar-d)/(mu_*(mu_*z+1.5*(zStar-d)))+2.0/3.0*log(1+1.5*(zStar-d)/(mu_*z))*z/(zStar-d));
    }

    inline const double* TAD_()
    {
	return &TAD;
    }

    inline double* LADTop_()
    {
	return &LADTop;
    }

    inline double uStar_()
    {
	return uStar;
    }

    inline double thetaStar_()
    {
	return thetaStar;
    }

    inline double L_()
    {
	return L;
    }

    inline double zc_()
    {
	return zc;
    }

    inline double zi_()
    {
	return zi;
    }

    inline double mu_()
    {
	return muM;
    }

    inline double Uh_(double uStar_, double L_)
    {
	return uStar_/(kappa*sparsityCorr)*(log((htree-d)/Z0) - stabKernel((htree-d)/L_,"M","integral")+stabKernel(Z0/L_,"M","integral"));
    }

    inline double GradUh_(double uStar_, double mu_, double L_)
    {
	return uStar_/(kappa*(htree-d))*(1-exp(-mu_*(htree-d)/(zStar-d)))*stabKernel((htree-d)/L_,"M");
    }

    inline double propConstTheory(double theta0)
    {
	double buoyancy = g/(theta0 + thetaStar/kappa*(log((zref-d)/zref_T0)-stabKernel((zref-d)/L,"H","integral")+stabKernel(zref_T0/L,"H","integral")+roughSublayerCorr(zref-d,muH,"H",L)-roughSublayerCorr(zref_T0,muH,"H",L)));
	return 1.4*pow(-uStar*thetaStar*buoyancy,1.0/3.0);
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
