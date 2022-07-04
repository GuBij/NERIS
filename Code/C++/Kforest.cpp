#include <iomanip>
#include <cmath>
#include "Kforest.h"

 Kforest::Kforest(const dictionary& dict)
 :
 forest(dict), 
 name_("Kforest"),
 ofilename("output/PARAMS_" + name_ + "_" + meteofile())
 {
 }

 void Kforest::K(double* k, const double& z)
 {
	double sigma_w_new;
        if (L < 0 && L > -200 && z > zStar)
          sigma_w_new = propConst*pow(z-forest::d_(),1.0/3.0);
        else if (L < 0 && L > -200 && z > htree)
          sigma_w_new = propConst*pow(zStar-forest::d_(),1.0/3.0);
        else if (z > htree)
          sigma_w_new = sigma_w;
        else if (L < 0 && L > -200)
        {
          sigma_w_new = propConst*pow(zStar-forest::d_(),1.0/3.0);
          sigma_w_new = 0.9*sigma_w_new/htree*z+0.1*sigma_w_new;
        }
        else if (z < htree)
          sigma_w_new = 0.9*sigma_w/htree*z+0.1*sigma_w;

        *k=0.5*sigma_u*sigma_u*timeStep_*(timeStep_*(2-timeStep_/tau_u)/tau_u);
        *(k+1)=0.5*sigma_v*sigma_v*timeStep_*(timeStep_*(2-timeStep_/tau_v)/tau_v);
        *(k+2)=0.5*sigma_w_new*sigma_w_new*vertTimeStep_*(vertTimeStep_*(2-vertTimeStep_/tau_w)/tau_w);
 }

 void Kforest::GradK(double* k, const double& z)
 {

	double Kmag;

	if (L < 0 && L > -200 && z > zStar)
	   Kmag = propConst*propConst*1.0/3.0*pow(1/(z-forest::d_()),1.0/3.0)*vertTimeStep_;
	else if (L < 0 && L > -200 && z < htree)
	{
	   Kmag = propConst*pow(zStar-forest::d_(),1.0/3.0);
	   Kmag = (0.9/htree*z+0.1)*Kmag*Kmag*0.9/htree*vertTimeStep_;
	}
	else if (z < htree)
	   Kmag = (0.9/htree*z+0.1)*sigma_w*sigma_w*0.9/htree*vertTimeStep_;

        *k=0.0;
        *(k+1)=0.0;
        *(k+2)=Kmag*(vertTimeStep_*(2-vertTimeStep_/tau_w)/tau_w);
 }

 void Kforest::updateParams(ifstream* ifile)
 {
	int timeIndex; 
        (*ifile) >> timeIndex >> sigma_u >> sigma_v >> sigma_w >> L >> propConst >> timeStep_ >> vertTimeStep_ >> tau_u >> tau_v >> tau_w;
 }

 void Kforest::writeProfile()
 {
        ofstream ofile("output/KList.txt", ofstream::trunc);

        double u[3], dz = 0.5;
        for (int i = 1; i < 200; i++)
        {
          K(u,i*dz);
          ofile << i*dz << setw(15) << u[0] << endl;
        }
        ofile.close();
 }
