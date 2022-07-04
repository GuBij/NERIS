#include <stdlib.h>
#include <algorithm>
#include <iterator>
#include "forest.h"
using namespace std;

 forest::forest(const dictionary& dict) 
 :
  dictionary(dict,"forest",{"LAI", "Cd", "htree", "canopyDepth", "zref", "standDensity", "DBH","Sct","zrefTemp","zrefTemp0"}),
  name_("forest"),
  plantCover(1-exp(-0.526*getfl(0))), 
  alpha(2*pow(getfl(1)*getfl(0),0.25)),
  htree(getfl(2)),
  zStar(2*getfl(2)),
  canopyDepth(getfl(3)),
  a2(4.1920),
  a3(1.9264),
  a4(1.0),
  a1(getfl(0)*(a2+1)/(gsl_sf_hyperg_2F1(a2+1,-a3,a2+2,1/a4)*getfl(3))),
  zref(getfl(4)),
  zref_T(getfl(8)),
  zref_T0(getfl(9)),
  kappa(0.4),
  RH(0.8),
  pressure(101300.0),
  rho(1.225),
  R(287.0),
  cp(1004.0),
  g(9.81),
  z0(0.071*getfl(2)),
  z0g(0.1),
  TAD(getfl(5)*0.5*acos(-1)*getfl(6)*1e-04),
  stackHeight(dict.getfl(2)),
  diam(dict.getfl(3)),
  flowRate(dict.getfl(4)),
  Ts(dict.getfl(5)+273.15),
  meteofile_(dict.get(0))
 {
   sparsityCorr = plantCover*(alpha-1)+1;
   d = plantCover*alpha*((a2*a4/(a2+a3)-1)*canopyDepth+htree)/sparsityCorr;
   Z0 = alpha*alpha*z0/sparsityCorr;

   ifstream ifile(meteofile_);
   ifile.unsetf(ios_base::skipws);
   if (ifile.is_open())
     noRecords = count(istream_iterator<char>(ifile),istream_iterator<char>(),'\n');
   else
     cerr << "file " << meteofile_ << " not found. Command executed from class forest." << endl;

   ifile.setf(ios_base::skipws);
   ifile.close();
 }

 double forest::stabKernel(double z, string type, string option)
 {
   if (option == "integral")
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
   else if (option == "derivative")
   {
      if (z<0 && (type == "M"))
	return 0.25*16*pow(1-16*z,-1.25);
      else if (z<0 && (type == "H"))
	return 16*0.5*pow(1-16*z,-1.5);
      else if (z>0)
	return 5.0;
      else
	return 0.0;
   }
 }

 double forest::stabKernel(double z, string type)
 {
     if (z > 0)
       return 1.0+5.0*z;
     else if (type == "H")
       return pow(1.0-16.0*z,-0.5);
     else if (type == "M")
       return pow(1.0-16.0*z,-0.25);
 }

 void forest::stateEqn(double* x, double* F)
 {
   F[0] = (log((htree-d)/z0) + roughSublayerCorr(htree-d,x[0]) - roughSublayerCorr(z0,x[0]) - log(sparsityCorr*(htree-d)/(alpha*alpha*z0))/sparsityCorr);
   F[0] *= -F[0];
 }

 void forest::JacStateEqn(double* x, double* J)
 {

   J[0] = 2*(log((htree-d)/z0) + roughSublayerCorr(htree-d,x[0]) - roughSublayerCorr(z0,x[0]) - log(sparsityCorr*(htree-d)/(alpha*alpha*z0))/sparsityCorr);
   J[0] *= JacRoughSublayerCorr(htree-d,x[0]) - JacRoughSublayerCorr(z0,x[0]); 
 };

 double forest::findMu(double mu0)
 {
    double f_SMALL_FLOAT = pow(2.0,-64.0);
    int MAX_NOIT = 10000;

    double fres=1, mu = mu0, lambda = 1.5;
    int noit = 0;

    double psiStar, psiStar0, J, F;
    while(abs(fres) > f_SMALL_FLOAT && noit < MAX_NOIT)
    {
	psiStar = roughSublayerCorr(htree-d,mu);
	psiStar0 = roughSublayerCorr(z0,mu);

	J = -(htree-d)/(zStar-d)*psiStar-(mu*(htree-d)/(mu*(htree-d)+lambda*(zStar-d)))*(zStar-d)/(mu*mu*(htree-d))*exp(-mu*(htree-d)/(zStar-d));
	J = J + z0/(zStar-d)*psiStar0+(mu*z0/(mu*z0+lambda*(zStar-d)))*(zStar-d)/(mu*mu*z0)*exp(-mu*z0/(zStar-d));
	F = log((htree-d)/z0)+psiStar-psiStar0-1/sparsityCorr*log((htree-d)/Z0);
	fres = F/mu;
	if ( mu > F/J )
	  mu -= F/J;
	else
	  noit = MAX_NOIT;

	noit += 1;
    }

    return mu;
 }

 double forest::findMu(double mu0, string type, double L)
 {
    double f_SMALL_FLOAT = pow(2.0,-64.0);
    int MAX_NOIT = 10000;

    double fres=1, mu = mu0, lambda = 1.5, nu = 0.5;
    int noit = 0;

    double J, J0, F, x1, x0, phi, phi0;
    while(abs(fres) > f_SMALL_FLOAT && noit < MAX_NOIT)
    {
	x1 = (1+nu*(zStar-d)/(mu*(htree-d)))*(htree-d)/L;
	x0 = (1+nu*(zStar-d)/(mu*z0))*z0/L;
	phi = stabKernel(x1,type);
	phi0 = stabKernel(x0,type);

	J = -nu*(zStar-d)/(mu*mu*L)*stabKernel(x1,type,"derivative")-(htree-d)/(zStar-d)*phi; J *= log(1+lambda*(zStar-d)/(mu*(htree-d)))/lambda;
	J -= (zStar-d)/(mu*mu*(htree-d))*phi*mu*(htree-d)/(mu*(htree-d)+lambda*(zStar-d)); J *= exp(-mu*(htree-d)/(zStar-d));
	J0 = -nu*(zStar-d)/(mu*mu*L)*stabKernel(x0,type,"derivative")-z0/(zStar-d)*phi0; J0 *= log(1+lambda*(zStar-d)/(mu*z0))/lambda;
	J0 -= (zStar-d)/(mu*mu*z0)*phi0*mu*z0/(mu*z0+lambda*(zStar-d)); J0 *= exp(-mu*z0/(zStar-d));
	J -= J0;
        F = log((htree-d)/z0)-stabKernel((htree-d)/L,type,"integral")+stabKernel(z0/L,type,"integral")+roughSublayerCorr(htree-d,mu,type,L)-roughSublayerCorr(z0,mu,type,L);
	F -= 1/sparsityCorr*(log((htree-d)/Z0)-stabKernel((htree-d)/L,type,"integral")+stabKernel(Z0/L,type,"integral"));
        fres = F/mu;
        if (mu > F/J)
          mu -= F/J;
        else
          noit = MAX_NOIT;

	noit += 1;
    }

    return mu;
 }

 double forest::picard(double Uref, double theta, double theta0, double* x0)
 {
   double f_SMALL_FLOAT = 1e-04;
   double MAX_NOIT = 10000;

   double muM = x0[0], uStar = x0[1], muH = x0[2], thetaStar = x0[3], L = x0[4], res = 1.0;
   double z1 = zref-d, z2 = zref_T - d, x[5];
   int noit = 0;
   while(res > f_SMALL_FLOAT && noit < MAX_NOIT)
   {
	muM=findMu(muM,"M",L);

	uStar = Uref/(log(z1/z0)-stabKernel(z1/L,"M","integral")+stabKernel(z0/L,"M","integral")+roughSublayerCorr(z1,muM,"M",L)-roughSublayerCorr(z0,muM,"M",L));

	muH=findMu(muH,"H",L);

	thetaStar = (theta-theta0)/(log(z2/z0)-stabKernel(z2/L,"H","integral")+stabKernel(z0/L,"H","integral")+roughSublayerCorr(z2,muH,"H",L)-roughSublayerCorr(z0,muH,"H",L));

	L = uStar*uStar*theta0/(g*thetaStar);

	x[0] = muM; x[1] = uStar; x[2] = muH; x[3] = thetaStar; x[4] = L;
	res = (x[0]-x0[0])*(x[0]-x0[0]); double min_x0 = x0[0]*x0[0];
	for (int i = 1; i < 5; i++)
	{
	  if ((x[i]-x0[i])*(x[i]-x0[i]) > res)
	    res = (x[i]-x0[i])*(x[i]-x0[i]);

	  if (x0[i]*x0[i] < min_x0)
	    min_x0 = x0[i]*x0[i];
	}
	res = sqrt(res/min_x0);
	noit += 1;
	x0[0] = muM; x0[1] = uStar; x0[2] = muH; x0[3] = thetaStar; x0[4] = L;
   }

   cout << "noit: " << noit << ", res: " << res << endl;
 }

 int forest::calibrate(double theta0, double theta, double Uref, function<void(double*,double*)> JacLAD_1, function<void(double*,double*)> JacLAD_2, function<void(double*,double*)> LAD_1, function<void(double*,double*)> LAD_2)
 {
   double x[1]; int status;
   x[0] = 0.25;
   status = math::findZero(x, 1, JacLAD_1, LAD_1);
   if (status == 0)
	zi = (x[0]-1)*canopyDepth+htree;
   else
	return status;

   muM = findMu(2.59);
   uStar = Uref/(log((zref-d)/z0)+roughSublayerCorr(zref-d,muM)-roughSublayerCorr(z0,muM)); 
   muH = muM;
   thetaStar = (theta-theta0)/(log((zref_T-d)/z0)+roughSublayerCorr(zref_T-d,muH)-roughSublayerCorr(z0,muH)); 
   L = uStar*uStar*theta0/(g*thetaStar); 
   double y[5] = {muM,uStar,muH,thetaStar,L};
   picard(Uref,theta,theta0,y);
   muM = y[0]; uStar = y[1]*kappa; muH = y[2]; thetaStar = y[3]*kappa; L = y[4];

   Uh = Uh_(uStar,L);
   GradUh = GradUh_(uStar, muM, L);
   LADTop = GradUh/(Uh*plantCover);
   x[0] = 0.95;
   status = math::findZero(x, 1, JacLAD_2, LAD_2);
   if (status == 0)
   {
	  zc = (x[0]-1)*canopyDepth+htree;
   }
   return status;
 }

 void forest::plumeRise(double U, double T, double T0, double uStar, double L, double& heff)
 {
   double dT = T - T0;
   double Tair = T0 + dT/(zref_T-zref_T0)*(stackHeight-zref_T0) + 273.15;
   double v = flowRate/(acos(-1)*0.5*0.5*diam*diam);
   double F = 0.25*g*v*diam*diam*(Ts-Tair)/Ts;
   double s = g*(dT/(zref_T-zref_T0)+0.0098)/Tair;

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

   heff += stackHeight;
 }

 void forest::LagrangianTimeScale(double heff, double L, double uStar, double thetaStar, double sigma_u, double sigma_v, double sigma_w, double* tau, double theta0, double& h)
 {
    double z;
    if (heff > htree)
	z = heff-d;
    else
	z = htree-d;

    if (abs(L)>200)
    {
        *tau = 0.5*z/(sigma_w*(1+30*7.2921e-5*sin(58*acos(-1)/180)*z/uStar));
        *(tau+1) = *tau; *(tau+2) = *tau;
    }
    else if (L > 0)
    {
        h = 0.25*sqrt(0.5*uStar*L/(7.2921e-5*sin(58*acos(-1)/180)))+d;
        *tau = 0.15*h/sigma_u*sqrt(z/h);
        *(tau+1) = 0.07*h/sigma_v*sqrt(z/h);
        *(tau+2) = 0.1*h/sigma_w*sqrt(z/h);
    }
    else
    {
	if (L > -100 || sigma_v < uStar*pow(12,1.0/3.0))
	{
	  double buoyancy = g/(theta0 + thetaStar/kappa*(log((zref-d)/zref_T0)-stabKernel((zref-d)/L,"H","integral")+stabKernel(zref_T0/L,"H","integral")+roughSublayerCorr(zref-d,muH,"H",L)-roughSublayerCorr(zref_T0,muH,"H",L))); 
          double wStar = sigma_v/0.6;
	  h = pow(wStar,3.0)/(-uStar*thetaStar*buoyancy)+d;	
	}
        else
          h = (pow(sigma_v/uStar,3.0)-12)*2*(-L);

        *(tau+1) = 0.15*h/sigma_v;
        *tau = *(tau+1);

        if (z/h > 0.1)
          *(tau+2) = 0.15*h/sigma_w*(1-exp(-5*z/h));
        else if (z/h < 0.1 && (z-z0) > -L)
          *(tau+2) = 0.1*z/(sigma_w*(0.55-0.38*(z-z0)/L));
        else
          *(tau+2) = 0.59*z/sigma_w;
    }
 }

 void forest::timeStepSelection(const int dim, const double* tauL, const double windSpeed, int* timeStep)
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
	   int m = 0, newTimeStep = 0; double min_dist = 1e+6; //bool notfound = true;
           while ( m < 21)
           {
		double dist = abs(candidates[m]-timeStep[2]);
                if ((timeStep[0] % candidates[m]) == 0 && candidates[m] > 0.02*tauL[2] && candidates[m] < 0.5*tauL[2] && dist < min_dist)
                {
                   //notfound = false;
                   newTimeStep = candidates[m];
		   min_dist = dist;
                }

                m++;
           }
          *(timeStep+2) = newTimeStep;
        }
   }
 }

 void forest::printParams()
 {
   cout << "\n\tmu: " << muM
	<< "\n\td: " << d
	<< "\n\tz0: " << z0
	<< "\n\tuStar: " << uStar 
	<< "\n\tzc: " << zc
	<< "\n\tzi: " << zi << "\n" << endl;
 }

