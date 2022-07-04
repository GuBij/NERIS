#include "KTaylorForest.h"

 KTaylorForest::KTaylorForest(const dictionary& dict)
 :
 forest(dict),
 name_("KTaylorForest")
 {
 }

 void KTaylorForest::K ( double* k, const double& z)
 {
	double sigw;
        if (L < 0 && L > -200 && z > zStar)
          sigw = propConst*pow(z-forest::d_(),1.0/3.0);
	else if (L < 0 && L > -200 && z > htree)
	  sigw = propConst*pow(zStar-forest::d_(),1.0/3.0);
	else if (z > htree)
	  sigw = sd_w;
	else if (L < 0 && L > -200)
	{
	  sigw = propConst*pow(zStar-forest::d_(),1.0/3.0);
	  sigw = 0.9*sigw/htree*z+0.1*sigw;
	}
	else
	  sigw = 0.9*sd_w/htree*z+0.1*sd_w;

        *k=sqrt(1-autoCorr[0]*autoCorr[0])*sd_u;
        *(k+1)=sqrt(1-autoCorr[1]*autoCorr[1])*sd_v;
        *(k+2)=sqrt(1-autoCorr[2]*autoCorr[2])*sigw;
 }

 double KTaylorForest::GradSig2 ( const double& z)
 {
	if (L < 0 && L > -200 && z > zStar)
	  return propConst*propConst*2.0/3.0*pow(1/(z-forest::d_()),1.0/3.0);
	else if (z > htree)
	  return 0.0;
	else
	  return 1.8*(0.9*sd_w/htree*z+0.1*sd_w)*sd_w/htree;
 }

 void KTaylorForest::updateParams (ifstream* ifile)
 {
        int timeIndex; double dummy, dummy2;
        (*ifile) >> timeIndex >> sd_u >> sd_v >> sd_w >> dummy >> tau_u >> tau_v >> tau_w >> propConst >> hABL >> L >> timeStep >> timeStepVert; 

	autoCorr[0] = 1-timeStep/tau_u;
	autoCorr[1] = 1-timeStep/tau_v;
	autoCorr[2] = 1-timeStepVert/tau_w;
 }
