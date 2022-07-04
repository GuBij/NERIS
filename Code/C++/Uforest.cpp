#include <iomanip>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>
#include "Uforest.h"

 Uforest::Uforest(const dictionary& dict)
 :
 Uforest(dict,false)
 {}

 Uforest::Uforest(const dictionary& dict, bool LangevinModel)
 :
 forest(dict), 
 name_("Uforest"),
 ofilename("output/PARAMS_" + name_ + "_" + meteofile()),
 Kfilename("output/PARAMS_Kforest_" + meteofile()),
 stackHeight(dict.getfl(2)),
 useLangevin(LangevinModel)
 {
   cout << "\nforest model selected: \n"
        << "\t - plant cover: " << plantCover << "\n"
        << "\t - tree height: " << htree << "\n"
        << "\t - thickness roughness layer: " << zStar << "\n"
        << "\t - thickness crown: " << canopyDepth << "\n"
        << "\t - normalization factor LAD model: " << a1 << "\n"
        << "\t - roughness length: " << z0 << "\n"
        << "\t - tree area density: " << TAD << "\n"
        << "\t - displacement height: " << forest::d_() << "\n" << endl;
 }

 void Uforest::calcParams(int timeIndex,double* data)
 {
   double T0 = data[0];
   double T = data[1];
   double uref=data[2];
   double windDirection_ = data[3];
   double elevation_ = data[4];

   const double eps=0.622;

   const double es0 = (1.0007+3.46*1e-06*pressure)*611.21*exp(17.502*T0/(240.97+T0)); 
   const double rv0 = eps*RH*es0/(pressure-RH*es0);
   const double theta0 = (T0+273.15)*(1+rv0/eps)/(1+rv0)*pow(pressure/(pressure-g*rho*zref_T0),R/cp);

   const double es = (1.0007+3.46*1e-06*pressure)*611.21*exp(17.502*T/(240.97+T)); 
   const double rv = eps*RH*es/(pressure-RH*es);
   const double theta = (T+273.15)*(1+rv/eps)/(1+rv)*pow(pressure/(pressure-g*rho*zref_T),R/cp);

   auto LAD2_1 = bind(&forest::LAD2,*this,placeholders::_1,placeholders::_2,TAD_());
   auto LAD2_2 = bind(&forest::LAD2,*this,placeholders::_1,placeholders::_2,LADTop_());
   auto JacLAD2_1 = bind(&forest::JacLAD2,*this,placeholders::_1,placeholders::_2,TAD_());
   auto JacLAD2_2 = bind(&forest::JacLAD2,*this,placeholders::_1,placeholders::_2,LADTop_());
   int status = calibrate(theta0,theta,uref,JacLAD2_1,JacLAD2_2,LAD2_1,LAD2_2);

   updateParams(forest::uStar_(), forest::mu_(), forest::zc_(), forest::zi_(), forest::L_(),windDirection_);
   double Ustack[3]; U(Ustack,stackHeight); Ustack[0] = sqrt(Ustack[0]*Ustack[0]+Ustack[1]*Ustack[1]+Ustack[2]*Ustack[2]);
   if (stackHeight > htree)
     forest::plumeRise(Ustack[0],T,T0,forest::uStar_(),L,heff_);
   else
     heff_ = stackHeight;
   int width=12;
   ofstream ofile(ofilename, ios::app);
   if (ofile.is_open())
   {
	ofile << left << setw(width) << timeIndex 
	     << left << setw(width) << forest::uStar_()
	     << left << setw(width) << forest::mu_()
	     << left << setw(width) << forest::zc_()
             << left << setw(width) << forest::zi_()
	     << left << setw(width) << windDirection_
	     << left << setw(width) << elevation_
	     << left << setw(width) << L
	     << left << setw(width) << heff_ << "\n";
	ofile.close();
   }

     double pi = acos(-1);
     double sigmaAz = data[5];
     double sigmaEl = data[6];
     double sigma_v = uref*sigmaAz*pi/180; 
     double sigma_w = uref*sigmaEl*pi/180; 

     double sigma_u;
     if (abs(L) > 200)
        sigma_u = 1.2*sigma_v;
     else if (L < 0 && L > -200)
        sigma_u = sigma_v;
     else
        sigma_u = 1.5*sigma_v;

     double tau[3], bound[2], hABL = 0.0, propConst = 0.0, propConstTheoretical = 0.0;
     forest::LagrangianTimeScale(heff_,L, forest::uStar_(), forest::thetaStar_(), sigma_u, sigma_v, sigma_w, tau, theta0, hABL);

     if (L<0 && L > -200)
     {
        propConst = sigma_w/pow((zref-forest::d_()),1.0/3.0);
	propConstTheoretical = forest::propConstTheory(theta0);
     }

     int timeStep[3];
     double Urelease[3]; U(Urelease,heff_); Urelease[0] = sqrt(Urelease[0]*Urelease[0]+Urelease[1]*Urelease[1]+Urelease[2]*Urelease[2]);
     forest::timeStepSelection(3,tau,Urelease[0],timeStep);

     if (timeStep[2] > timeStep[0])
     {
        cerr << "WARNING: vertical time step is bigger than horizontal time step in " << meteofile() << "\n" << endl;
	exit(EXIT_FAILURE);
     }
     else if (timeStep[0] <= 0 || timeStep[2] <= 0)
     {
        cerr << "WARNING: time step is zero in " << meteofile() << "\n" << endl;
	exit(EXIT_FAILURE);
     }

   if (useLangevin)
   {
     ofstream Kfile("output/PARAMS_KTaylorForest_" + meteofile(), ios::app);
     if (Kfile.is_open())
     {
        Kfile << left << setw(width) << timeIndex
              << left << setw(width) << sigma_u
	      << left << setw(width) << sigma_v
              << left << setw(width) << sigma_w
              << left << setw(width) << forest::uStar_()
              << left << setw(width) << tau[0]
              << left << setw(width) << tau[1]
              << left << setw(width) << tau[2] 
              << left << setw(width) << propConst
              << left << setw(width) << hABL
              << left << setw(width) << L 
	      << left << setw(width) << timeStep[0]
              << left << setw(width) << timeStep[2] << "\n";
        Kfile.close();
     }
   }
   else
   {
     ofstream Kfile(Kfilename, ios::app);
     if (Kfile.is_open())
     {
        int width=15;
        Kfile << left << setw(width) << timeIndex 
              << left << setw(width) << sigma_u
              << left << setw(width) << sigma_v
              << left << setw(width) << sigma_w
              << left << setw(width) << L
              << left << setw(width) << propConst
              << left << setw(width) << timeStep[0]
              << left << setw(width) << timeStep[2] 
              << left << setw(width) << tau[0]
              << left << setw(width) << tau[1]
              << left << setw(width) << tau[2] << "\n";
        Kfile.close();
     }
   }
 }

 double Uforest::intFGF(double z)
 {
	double intFGF_value;

	if (z > htree)
	  intFGF_value = 0.0;
	else if (z > zc)
	  intFGF_value = (htree-z)*GradUh/(plantCover*Uh);
	else if (z > zi)
	{
	  double y=(z-htree)/canopyDepth+1;
	  intFGF_value = dummy[1]-a1*canopyDepth/(a2+1)*pow(y,a2+1)*pow(a4,a3)*gsl_sf_hyperg_2F1(a2+1,-a3,a2+2,y/a4);
	}
	else if (z > 10*z0g)
	  intFGF_value = dummy[2]+TAD*(zi-z);
	else
 	{
	  intFGF_value = dummy[3]+10*z0g*TAD*log(10*z0g/z);
	}

	return intFGF_value;
 }

 void Uforest::U(double* u, const double& z)
 {
	double Umag;

        if (z > htree)
          Umag = uStar/kappa*(log((z-forest::d_())/z0) - forest::stabKernel((z-forest::d_())/L,"M","integral") + forest::stabKernel(z0/L,"M","integral") + roughSublayerCorr(z-forest::d_(),mu,"M",L) + dummy[0]); 
        else
          Umag = Uh*exp(-plantCover*intFGF(z));

        *u=Umag*sin(windDirection); 
        *(u+1)=Umag*cos(windDirection); 
        *(u+2)=0.0;
 }

 void Uforest::updateParams(ifstream* ifile) 
 {
	int timeIndex; 
	  *ifile >> timeIndex >> uStar >> mu >> zc >> zi >> windDirection >> elevation >> L >> heff_;

	  Uh = Uh_(uStar,L);
	  GradUh = GradUh_(uStar,mu,L);
	  elevation *= pi/180; windDirection *= pi/180; windDirection += pi;

          double yc=(zc-htree)/canopyDepth+1, F21c = pow(yc,a2+1)*pow(a4,a3)*gsl_sf_hyperg_2F1(a2+1,-a3,a2+2,yc/a4);
          double yi=(zi-htree)/canopyDepth+1, F21i = pow(yi,a2+1)*pow(a4,a3)*gsl_sf_hyperg_2F1(a2+1,-a3,a2+2,yi/a4);
	  dummy[0]=-roughSublayerCorr(z0,mu,"M",L);
	  dummy[1]=(htree-zc)*GradUh/(plantCover*Uh)+a1*canopyDepth/(a2+1)*F21c;
	  dummy[2]=dummy[1]-a1*canopyDepth/(a2+1)*F21i;
	  dummy[3]=dummy[2]+TAD*(zi-10*z0g);
 }

 void Uforest::updateParams(double uStar_, double mu_, double zc_, double zi_, double L_, double windDirection_) 
 {
	  uStar = uStar_; mu = mu_; zc = zc_; zi = zi_; windDirection = windDirection_; L = L_;

          Uh = Uh_(uStar,L);
          GradUh = GradUh_(uStar,mu,L);
          elevation *= pi/180; windDirection *= pi/180; windDirection += pi;

          double yc=(zc-htree)/canopyDepth+1, F21c = pow(yc,a2+1)*pow(a4,a3)*gsl_sf_hyperg_2F1(a2+1,-a3,a2+2,yc/a4);
          double yi=(zi-htree)/canopyDepth+1, F21i = pow(yi,a2+1)*pow(a4,a3)*gsl_sf_hyperg_2F1(a2+1,-a3,a2+2,yi/a4);
          dummy[0]=-roughSublayerCorr(z0,mu,"M",L);
          dummy[1]=(htree-zc)*GradUh/(plantCover*Uh)+a1*canopyDepth/(a2+1)*F21c;
          dummy[2]=dummy[1]-a1*canopyDepth/(a2+1)*F21i;
          dummy[3]=dummy[2]+TAD*(zi-10*z0g);
 }

 void Uforest::writeProfile()
 {
	ofstream ofile("output/UProfile.txt", ofstream::trunc);

	double u[3], dz = 0.5;
	for (int i = 1; i < 200; i++)
	{
	  U(u,i*dz);
	  ofile << i*dz << setw(15) << u[0] << endl;
	}
	ofile.close();
 }
