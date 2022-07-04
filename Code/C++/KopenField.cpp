#include "KopenField.h"

 KopenField::KopenField(const dictionary& dict)
 :
 openField(dict),
 name_("KopenField")
 {
 }

 void KopenField::K ( double* k, const double& z)
 {
        if (L < 0 && L > -200)
          sigma_w = propConst*pow(z,1.0/3.0);

        *k=0.5*sigma_u*sigma_u*timeStep_*(timeStep_*(2-timeStep_/tau_u)/tau_u);
        *(k+1)=0.5*sigma_v*sigma_v*timeStep_*(timeStep_*(2-timeStep_/tau_v)/tau_v);
        *(k+2)=0.5*sigma_w*sigma_w*vertTimeStep_*(vertTimeStep_*(2-vertTimeStep_/tau_w)/tau_w);
 }

 void KopenField::GradK ( double* k, const double& z)
 {
	*k=0;
	*(k+1)=0;

        if (L < 0 && L > -200 && z < 2*z0)
          *(k+2) = propConst*propConst*1.0/3.0*pow(0.5/z0,1.0/3.0)*vertTimeStep_*(vertTimeStep_*(2-vertTimeStep_/tau_w)/tau_w);
        else if (L< 0 && L > -200)
          *(k+2) = propConst*propConst*1.0/3.0*pow(1/z,1.0/3.0)*vertTimeStep_*(vertTimeStep_*(2-vertTimeStep_/tau_w)/tau_w);
	else
	  *(k+2) = 0;
 }

 void KopenField::updateParams (ifstream* ifile)
 {
        int timeIndex;
        (*ifile) >> timeIndex >> sigma_u >> sigma_v >> sigma_w >> L >> propConst >> timeStep_ >> vertTimeStep_ >> tau_u >> tau_v >> tau_w;
 }
