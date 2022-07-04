 #include <iostream>
 #include <iomanip>
 #include <cmath>
 #include <fstream>
 #include <sstream>
 #include "radiation.h"

 radiation::radiation(int noRecords_)
 :
 E(1.294),
 A(59.4),
 L(12.0),
 mu(0),
 a1(0),
 a2(0),
 sigair(0),
 noRecords(noRecords_)
 {
   ifstream msfile("monitoringStations.txt");
   if (msfile.is_open())
   {
     double x, y;
     bool no_eof=true;
     while (no_eof)
     {
	msfile >> x >> y;
	if (!msfile.eof())
	  location.push_back({x,y});
	else
	  no_eof=false;
     }
     msfile.close();
   }
   else
   {
     location={{200.143503293714,212.237160795659},{200.155246861131,211.962463226418},{199.982903434886,211.936320720149},{199.855503938484,211.93776306347},{199.869570903673,212.167760330885},{199.904281147205,212.263778117704},{199.961095905096,212.29068811741}};

     for ( int i=0; i < 7; i++)
     {
        get<0>(location[i])=(get<0>(location[i])-199.986)*1000; get<1>(location[i])=(get<1>(location[i])-212.124)*1000.0;
     }
   }
   nor = location.size();
   cout << "\n Dose rates will be calculated for " << nor << " monitoring stations.\n" << endl;

   weightCoeff = new double[noRecords*noRecords*nor]();
   neighbWeightCoeff = new double[noRecords*noRecords*nor]();
   dose = new double[nor*noRecords]();
   neighb = new double[nor*noRecords]();
 
  double logE=log(E);
  vector<double> A1={-8.46767*pow(10,-5),2.71806*pow(10,-4),2.93831*pow(10,-3),-6.30793*pow(10,-3),-1.84327*pow(10,-2),4.78301*pow(10,-2),-4.61222*pow(10,-2)};
  vector<double> A2={-7.15631*pow(10,-4),-2.54597*pow(10,-3),9.31948*pow(10,-3),1.89129*pow(10,-2),-3.12917*pow(10,-2),1.49122*pow(10,-2),-1.02463*pow(10,-2)};
  vector<double> A3={-1.40465*pow(10,-4),3.07113*pow(10,-4),1.42298*pow(10,-2),-3.57795*pow(10,-3),-1.18921*pow(10,-1),-4.20821*pow(10,-1),-2.70365};

  for (int i=0; i<7; i++)
  {
    a1+=A1[i]*pow(logE,6-i); a2+=A2[i]*pow(logE,6-i); mu+=A3[i]*pow(logE,6-i);
  }
  mu=exp(mu)*0.12041;

  if (E<0.1)
  {
    vector<double> a={-9.27823*pow(10,-3),0.262726,2.92391,15.8909,38.9774,30.6711};
    for ( int i=0 ;i < 6; i++)
	sigair+=a[i]*pow(logE,5-i);

    sigair=exp(sigair)/10.0;
  }
  else
  {
    vector<double> a={1.06134*pow(10,-5),-1.59201*0.001,7.51537*0.001,2.35368*0.01,-0.119158,-0.182464,-3.57954};
    for ( int i=0; i < 7; i++)
	sigair+=a[i]*pow(logE,6-i);

    sigair=exp(sigair)/10.0;
  }
 }

 void radiation::update(double* pos, int time, double Qp)
 {
   double rjsqr, mj;
   for ( int j = 0; j < nor; j++)
   {
	rjsqr = pow(pos[0]-get<0>(location[j]),2)+pow(pos[1]-get<1>(location[j]),2)+pow(pos[2]-1.5,2); 
	if ( rjsqr < (L*L))
	{
	  neighb[nor*(time - 1) + j] += Qp; 
	}
	else
	{
	  mj=mu*sqrt(rjsqr); dose[nor*(time - 1) + j] += Qp*(A*exp(-a1*mj)+(1-A)*exp(-a2*mj))/rjsqr*exp(-mj); 
	}
   }
 }

 void radiation::updateBayOpt(double* pos, int time, double Qp, int release, double scale)
 {
   double rjsqr, mj;
   for ( int j = 0; j < nor; j++)
   {
        rjsqr = pow(pos[0]-get<0>(location[j]),2)+pow(pos[1]-get<1>(location[j]),2)+pow(pos[2]-1.5,2);
        if ( rjsqr < (L*L))
        {
          neighb[nor*(time - 1) + j] += Qp; 
	  neighbWeightCoeff[(time-1)*noRecords*nor+nor*release+j] += scale;
        }
        else
        {
          mj=mu*sqrt(rjsqr); double doseKernel = (A*exp(-a1*mj)+(1-A)*exp(-a2*mj))/rjsqr*exp(-mj);
          dose[nor*(time - 1) + j] += Qp*doseKernel; 

          weightCoeff[(time-1)*noRecords*nor+nor*release+j] += scale*doseKernel;
        }
   }
 }

 void radiation::write(string meteofile)
 {
   double V = acos(-1)*(2.0*L*L*L/3.0 + L*L*1.5-1.5*1.5*1.5/3.0);
   double K=1.6*pow(10,-13+9), A1=(a1+1)*mu, A2=(a2+1)*mu, doseint = (-A*exp(-A1*L)/A1 - (1-A)*exp(-A2*L)/A2 + A/A1 + (1-A)/A2)*0.5*K*sigair*E/V; 

   ofstream ofile("output/Dose_" + meteofile); ofile.precision(6);
   for ( int i = 0; i < noRecords; i++)
   {
	int k=nor*i;
	for ( int j = 0; j < nor; j++)
	{
	  dose[k+j] *= 0.25*K*sigair*E/acos(-1); dose[k+j] += neighb[k+j]*doseint;
	  ofile << left << setw(15) << dose[k+j];
	}
	ofile << endl;
   }

   ofile.close();
 }

 void radiation::writeCoeffBayOpt(string meteofile)
 {
   double V = acos(-1)*(2.0*L*L*L/3.0 + L*L*1.5-1.5*1.5*1.5/3.0);
   double K=1.6*pow(10,-13+9), A1=(a1+1)*mu, A2=(a2+1)*mu, doseint = (-A*exp(-A1*L)/A1 - (1-A)*exp(-A2*L)/A2 + A/A1 + (1-A)/A2)*0.5*K*sigair*E/V; 

   ofstream ofile("output/weights_" + meteofile); ofile.precision(6);
   for ( int l = 0; l < noRecords; l++)
   {
     int m = noRecords*nor*l;
     for ( int i = 0; i < noRecords; i++)
     {
        int k=nor*i;
        for ( int j = 0; j < nor; j++)
        {
          weightCoeff[k+j+m] *= 0.25*K*sigair*E/acos(-1); weightCoeff[k+j+m] += neighbWeightCoeff[k+j+m]*doseint;
          ofile << left << i << " \t" << setw(15) << j << setw(15) << l << setw(15) << weightCoeff[k+j+m] << endl;
        }
     }
   }

   ofile.close();
 }

 void radiation::reconstruct(int noProcs, string meteofile)
 {
   const int nor_=7, noRecords_=25;
   stringstream ss; string id; double reconstrDose[nor_*noRecords_]={};
   for ( int j = 0; j < noProcs; j++)
   {
     ss << j; ss >> id;
     ifstream ifile("processor" + id + "/output/Dose_" + meteofile);
     for ( int r=0; r < (nor_*noRecords_); r++)
     {
	double dose; ifile >> dose;
	reconstrDose[r] += dose;
     }
   }

   ofstream ofile("output/Dose_" + meteofile, ofstream::trunc); ofile.precision(6);
   ofile << "period\tIMR/M11\tIMR/M12\tIMR/M07\tIMR/M08\tIMR/M13\tIMR/M09\tIMR/M10" << endl;
   for ( int i = 0; i < noRecords_; i++)
   {
	ofile << i+1; 
        int k=nor_*i;
        for ( int j = 0; j < nor_; j++)
          ofile << "\t" << left << setw(15) << reconstrDose[k+j];

        ofile << endl;
   }
   ofile.close();
 }

 void radiation::reconstructWeightsBayOpt(int noProcs, string meteofile)
 {
   const int nor_=7, noRecords_=25;
   stringstream ss; string id; double reconstrWeights[nor_*noRecords_*noRecords_]={};
   for ( int j = 0; j < noProcs; j++)
   {
     ss << j; ss >> id;
     ifstream ifile("processor" + id + "/output/weights_" + meteofile);
     for ( int r=0; r < (nor_*noRecords_*noRecords_); r++)
     {
        int i,j,l; double weight; ifile >> i >> j >> l >> weight;
        reconstrWeights[r] += weight;
     }
   }

   ofstream ofile("output/weights_" + meteofile, ofstream::trunc); ofile.precision(6);
   for ( int l = 0; l < noRecords_; l++)
   {
     int m = noRecords_*nor_*l;
     for ( int i = 0; i < noRecords_; i++)
     {
        int k=nor_*i;
        for ( int j = 0; j < nor_; j++)
          ofile << "\t" << left << i << " \t" << setw(15) << j << setw(15) << l << setw(15) << reconstrWeights[k+j+m] << endl;
     }
   }
   ofile.close();
 }

