#include "KTaylorOF.h"

 KTaylorOF::KTaylorOF(const dictionary& dict)
 :
 openField(dict,true),
 name_("KTaylorOF")
 {
 }

 void KTaylorOF::K ( double* k, const double& z)
 {
	if (L < 0 && L > -200)
	  sd_w = propConst*pow(z,1.0/3.0);

        *k=sqrt(1-autoCorr[0]*autoCorr[0])*sd_u;
        *(k+1)=sqrt(1-autoCorr[1]*autoCorr[1])*sd_v;
        *(k+2)=sqrt(1-autoCorr[2]*autoCorr[2])*sd_w;
 }

  double KTaylorOF::GradSig2 ( const double& z)
 {
	if (L < 0 && L > -200 && z < 2*z0)
	  return propConst*propConst*2.0/3.0*pow(0.5/z0,1.0/3.0);
	else if (L< 0 && L > -200)
	  return propConst*propConst*2.0/3.0*pow(1/z,1.0/3.0);
	else
	  return 0.0;
 }

 void KTaylorOF::updateParams (ifstream* ifile)
 {
        int timeIndex; double dummy, dummy2;
        (*ifile) >> timeIndex >> sd_u >> sd_v >> sd_w >> dummy >> tau_u >> tau_v >> tau_w >> propConst >> hABL >> L >> timeStep >> timeStepVert;

	autoCorr[0] = 1-timeStep/tau_u;
	autoCorr[1] = 1-timeStep/tau_v;
	autoCorr[2] = 1-timeStepVert/tau_w;
 }
