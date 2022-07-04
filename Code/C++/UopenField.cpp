#include "UopenField.h"

 UopenField::UopenField(const dictionary& dict)
 :
 openField(dict),
 name_("UopenField")
 {}

 void UopenField::U(double* u, const double& z)
 {
        double umag;

	if (z < 2*z0)
	  umag = uStar/kappa*(log(2)-openField::stabKernel(2*z0/L,"M","integral")+openField::stabKernel(z0/L,"M","integral"))/(2*z0)*z;
	else
	  umag = uStar/kappa*(log(z/z0)-openField::stabKernel(z/L,"M","integral")+openField::stabKernel(z0/L,"M","integral"));

        *u=umag*sin(windDirection); 
        *(u+1)=umag*cos(windDirection); 
        *(u+2)=0.0;
 }

 void UopenField::updateParams(ifstream* ifile)
 {
	int timeIndex;
	*ifile >> timeIndex >> uStar >> L >> heff_ >> windDirection >> elevation;

	windDirection *= pi/180; windDirection += pi; elevation *= pi/180;
 }

