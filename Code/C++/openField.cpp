#include <iomanip>
#include <algorithm>
#include <iterator>
#include "openField.h"

 openField::openField(const dictionary& dict)
 :
 openField(dict,false)
 {}
 
 openField::openField(const dictionary& dict, bool LangevinModel)
 :
 dictionary(dict,"openField",{"z0","zref","Sct","zrefTemp","zrefTemp0"}),
 name_("openField"),
 kappa(0.4),
 g(9.81),
 z0(getfl(0)),
 zref(getfl(1)),
 zref_T(getfl(3)),
 zref_T0(getfl(4)),
 Sct(getfl(2)),
 RH(0.8),
 pressure(101300.0),
 rho(1.225),
 R(287.0),
 cp(1004.0),
 meteofile_(dict.get(0)),
 stackHeight(dict.getfl(2)),
 d(dict.getfl(3)),
 flowRate(dict.getfl(4)),
 Ts(dict.getfl(5)+273.15),
 useLangevin(LangevinModel)
 {
   ifstream ifile(dict.get(0));
   ifile.unsetf(ios_base::skipws);
   if (ifile.is_open())
     noRecords = count(istream_iterator<char>(ifile),istream_iterator<char>(),'\n');
   else
     cerr << "file " << dict.get(0) << "not found. Command executed from class openField." << endl;

   ifile.setf(ios_base::skipws);
   ifile.close();

   cout << "\nopen field model selected." << endl;
 }

 void openField::picard(double theta0, double theta, double uref,  double& uStar, double& thetaStar, double& L)
 {
    double x0[3], x[3];
    x0[0] = uref/log(zref/z0);
    x0[1] = (theta-theta0)/log(zref_T/zref_T0);
    x0[2] = x0[0]*x0[0]*theta0/(g*x[1]);

    int noit=0, MAX_NOIT = 1000;
    double res=1.0, f_SMALL_FLOAT = 1e-04;
    while ( res > f_SMALL_FLOAT && noit < MAX_NOIT)
    {
      x[0] = uref/(log(zref/z0)-stabKernel(zref/x0[2],"M",true)+stabKernel(z0/x0[2],"M",true));
      x[1] = (theta-theta0)/(log(zref_T/zref_T0)-stabKernel(zref_T/x0[2],"H",true)+stabKernel(zref_T0/x0[2],"H",true));
      x[2] = x[0]*x[0]*theta0/(g*x[1]);
      res = (x[0]-x0[0])*(x[0]-x0[0])+(x[1]-x0[1])*(x[1]-x0[1])+(x[2]-x0[2])*(x[2]-x0[2]);
      res = sqrt(res);
      x0[0] = x[0]; x0[1] = x[1]; x0[2] = x[2];
      noit += 1;
    }

    uStar = x[0]*kappa; thetaStar = x[1]*kappa; L = x[2];
 }

 void openField::plumeRise(double T, double T0, double uStar, double L, double& heff)
 {
   double dT = T - T0;
   double Tair = T0 + dT/(zref_T-zref_T0)*(stackHeight-zref_T0) + 273.15;
   double v = flowRate/(acos(-1)*0.5*0.5*d*d);
   double F = 0.25*g*v*d*d*(Ts-Tair)/Ts;
   double s = g*(dT/(zref_T-zref_T0)+0.0098)/Tair;
   double U = uStar/kappa*(log(stackHeight/z0)-stabKernel(stackHeight/L,"M",true)+stabKernel(z0/L,"M",true));

   double xf;
   if (F < 55)
    xf = 49*pow(pow(F,10.0),1.0/16.0);
   else
    xf = 119*pow(F,0.4);

   if (L < 500 && L > 0)
   {
     if (xf > 1.84*U/pow(s,0.5))
	heff = 2.4*pow(F/(U*s),1.0/3.0);
     else
	heff = 1.6*pow(F,1.0/3.0)*pow(xf,2.0/3.0)/U;
   }
   else
   {
     int sign; (F > 0) ? (sign=1) : ((F < 0) ? (sign=-1) : (sign=0));
     heff = 1.6*sign*pow(abs(F),1.0/3.0)*pow(xf,2.0/3.0)/U;
   }

   heff = heff+stackHeight;
 }

 void openField::LagrangianTimeScale(double heff, double L, double uStar, double thetaStar, double sigma_u, double sigma_v, double sigma_w, double* tau, double theta0, double& h)
 {
    if (abs(L)>200)
    {
	*tau = 0.5*heff/(sigma_w*(1+30*7.2921e-5*sin(58*acos(-1)/180)*heff/uStar));
	*(tau+1) = *tau; *(tau+2) = *tau;
    }
    else if (L > 0)
    {
	h = 0.25*sqrt(0.5*uStar*L/(7.2921e-5*sin(58*acos(-1)/180)));
	*tau = 0.15*h/sigma_u*sqrt(heff/h);
	*(tau+1) = 0.07*h/sigma_v*sqrt(heff/h);
	*(tau+2) = 0.1*h/sigma_w*sqrt(heff/h);
    }
    else
    {
	if (L > -100 || sigma_v < uStar*pow(12,1.0/3.0) || (pow(sigma_v/uStar,3.0)-12)*2*(-L) < 100)
	{
          double buoyancy = g/(theta0 + thetaStar/kappa*(log(zref/zref_T0)-stabKernel(zref/L,"H",true)+stabKernel(zref_T0/L,"H",true)));
	  double wStar = sigma_v/0.6;
	  h = pow(wStar,3.0)/(-uStar*thetaStar*buoyancy);
	}
	else
	  h = (pow(sigma_v/uStar,3.0)-12)*2*(-L);

	*(tau+1) = 0.15*h/sigma_v;
	*tau = *(tau+1);

	if (heff/h > 0.1)
	  *(tau+2) = 0.15*h/sigma_w*(1-exp(-5*heff/h));
	else if (heff/h < 0.1 && (heff-z0) > -L)
	  *(tau+2) = 0.1*heff/(sigma_w*(0.55-0.38*(heff-z0)/L));
	else
	  *(tau+2) = 0.59*heff/sigma_w;
    }
 }

 void openField::timeStepSelection(const int dim, const double* tauL, const double windSpeed, int* timeStep)
 {
   int candidates[21] = {1,2,3,4,6,8,10,12,15,20,24,25,30,40,50,60,75,100,150,200,300}; //divisors of 600

   int i = 1; double max_tauL = *tauL, min_tauL = *tauL;
   while ( i < dim )
   {
      if (max_tauL < *(tauL+i))
	max_tauL = *(tauL+i);

      if (min_tauL > *(tauL+i))
	min_tauL = *(tauL+i);

      i++;
   }

   i = 0; bool noCandidate = true;
   while ( (max_tauL > min_tauL) && noCandidate && i < 21)
   {
	if (candidates[i] > 0.02*max_tauL && candidates[i] < 0.5*min_tauL)
		noCandidate = false;

	i++;
   }

   if (min_tauL == max_tauL)
	noCandidate = false;

   if (0.02*max_tauL < 0.5*min_tauL && !noCandidate)
   {
      int n = 10;
      if (windSpeed > 1.5)
      {
	double minTimeLapse = 120/windSpeed;
        while ( (n*0.01*max_tauL > min_tauL/n || minTimeLapse < (20*n*0.01*max_tauL)) && n > 2 )
        {
          n--;
        }
      }
      else if (windSpeed > 0)
      {
        double minTimeLapse = 120/windSpeed;
        while ( (n*0.01*max_tauL > min_tauL/n || minTimeLapse < (40*n*0.01*max_tauL)) && n > 2 )
        {
          n--;
        }
      }
      else
      {
        while ( n*0.01*max_tauL > min_tauL/n && n > 2 )
        {
          n--;
        }
      }

      int k = 0; bool notfound = true; double min_dist = 1e+6;
      while ( k < 21 && notfound)
      {
	if (candidates[k] > n*0.01*max_tauL && candidates[k] < min_tauL/n)
	{
	   for (int j = 0; j < dim; j++)
	     *(timeStep+j) = candidates[k];

	   notfound = false;
	}
	else if ( candidates[k] > 0.02*max_tauL && candidates[k] < 0.5*min_tauL && abs(candidates[k]-n*0.01*max_tauL) < min_dist)
	{
           for (int j = 0; j < dim; j++)
             *(timeStep+j) = candidates[k];

	   min_dist =  abs(candidates[k]-n*0.01*max_tauL);
	}  

	k++;
      }
   }
   else
   {
	timeStepSelection(2,tauL,windSpeed,timeStep);
	timeStepSelection(1,tauL+2,0,timeStep+2);
	if ((timeStep[0] % timeStep[2]) > 0)
	{
	   int m = 0; int newTimeStep = 0; double min_dist = 1e+6;
	   while ( m < 21 )
	   {
		double dist = abs(candidates[m]-timeStep[2]);
		if ((timeStep[0] % candidates[m]) == 0 && candidates[m] > 0.02*tauL[2] && candidates[m] < 0.5*tauL[2] && dist < min_dist)
		{
		   newTimeStep = candidates[m];
		   min_dist = dist;
		}

		m++;
	   }

	  *(timeStep+2) = newTimeStep;
	}
   }
 }

 void openField::calcParams(int timeIndex,double* data)
 {
   double T0 = data[0];
   double T = data[1];
   double uref=data[2];
   double windDirection = data[3];
   double elevation = data[4];

   double eps=0.622;

   double es0 = (1.0007+3.46*1e-06*pressure)*611.21*exp(17.502*T0/(240.97+T0)); 
   double rv0 = eps*RH*es0/(pressure-RH*es0);
   double theta0 = (T0+273.15)*(1+rv0/eps)/(1+rv0)*pow(pressure/(pressure-g*rho*zref_T0),R/cp);

   double es = (1.0007+3.46*1e-06*pressure)*611.21*exp(17.502*T/(240.97+T)); 
   double rv = eps*RH*es/(pressure-RH*es);
   double theta = (T+273.15)*(1+rv/eps)/(1+rv)*pow(pressure/(pressure-g*rho*zref_T),R/cp);

   double uStar, thetaStar, L, heff, hABL = 0.0;
   picard(theta0, theta, uref, uStar, thetaStar, L);
   plumeRise(T,T0,uStar,L,heff);

   int width=10;
   ofstream Ufile("output/PARAMS_UopenField_" + meteofile_, ios::app);
   if (Ufile.is_open())
   {
        Ufile << left << setw(width) << timeIndex
             << left << setw(width) << uStar
	     << left << setw(width) << L 
	     << left << setw(width) << heff
             << left << setw(width) << windDirection
	     << left << setw(width) << elevation << "\n";
        Ufile.close();
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

     double tau[3], propConst = 0.0, propConstTheory = 0.0, bound[2];
     LagrangianTimeScale(heff, L, uStar, thetaStar, sigma_u, sigma_v, sigma_w, tau, theta0, hABL);

     if (L<0 && L > -200)
     {
	propConst = sigma_w/pow(zref,1.0/3.0);
        double buoyancy = g/(theta0 + thetaStar/kappa*(log(zref/zref_T0)-stabKernel(zref/L,"H","integral")+stabKernel(zref_T0/L,"H","integral")));
        propConstTheory = 1.4*pow(-uStar*thetaStar*buoyancy,1.0/3.0);
     }

     int timeStep[3];
     double Urelease = uStar/kappa*(log(heff/z0)-stabKernel(heff/L,"M",true)+stabKernel(z0/L,"M",true));
     timeStepSelection(3,tau,Urelease,timeStep);

     if (timeStep[2] > timeStep[0])
     {
	cerr << "ERROR: vertical time step is bigger than horizontal time step in " << meteofile_ << "\n" << endl;
	exit(EXIT_FAILURE);
     }
     else if (timeStep[0] <= 0 || timeStep[2] <= 0)
     {
	cerr << "ERROR: time step is zero in " << meteofile_ << "\n" << endl;
	exit(EXIT_FAILURE);
     }

   if (useLangevin)
   {
     ofstream Kfile("output/PARAMS_KTaylorOF_" + meteofile_, ios::app);
     if (Kfile.is_open())
     {
        Kfile << left << setw(width) << timeIndex
              << left << setw(width) << sigma_u
	      << left << setw(width) << sigma_v 
	      << left << setw(width) << sigma_w
	      << left << setw(width) << uStar 
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
     ofstream Kfile("output/PARAMS_KopenField_" + meteofile_, ios::app);
     if (Kfile.is_open())
     {
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

 double openField::stabKernel(double z, string type, bool integrate)
 {
   if (integrate)
   {
     if (z > 0)
       return -5.0*z;
     else if (type == "H")
     {
       double x = pow(1-16.0*z,0.25);
       return 2*log(0.5*(1+x*x));
     }
     else if (type == "M")
     {
       double x = pow(1-16.0*z,0.25);
       return log(0.5*(1+x*x)*0.25*(1+2*x+x*x))-2*atan(x)+0.5*acos(-1);
     }
   }
   else
   {
     if (z > 0)
       return 1.0+5.0*z;
     else if (type == "H")
       return pow(1.0-16.0*z,-0.5);
     else if (type == "M")
       return pow(1.0-16.0*z,-0.25);
   }
 }
