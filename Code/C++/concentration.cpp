 #include <iomanip>
 #include <fstream>
 #include <sstream>
 #include <cmath>
 #include "concentration.h"

 concentration::concentration(const dictionary& dict, meteorology* wind) 
 :
 dictionary(dict,"sampleConc",{"Xmin","Xmax","Ymin","Ymax","Zmin","Zmax","resX","resY","resZ"}),
 varPos(0),
 pfileName("output/Positions_" + wind->meteofile()),
 ofileName("output/CONC_" + wind->meteofile()),
 mfreq(dict.getfl(2)),
 ppr(ceil(dict.getfl(3)/dict.getfl(8))),
 xMax(dict.getfl(4)),
 yMax(dict.getfl(5)),
 zMax(dict.getfl(6)),
 Xmin(getfl(0)),
 Xmax(getfl(1)),
 Ymin(getfl(2)),
 Ymax(getfl(3)),
 Zmin(getfl(4)),
 Zmax(getfl(5)),
 resX(getfl(6)),
 resY(getfl(7)),
 resZ(getfl(8)),
 nop(0)
 {
   bandWidth.reserve(2000);
   position = new double[3*ppr];
   meanDrift[0] = 0; meanDrift[1] = 0; meanDrift[2] = 0;
   ifstream windParamFile("output/PARAMS_" + wind->name() + "_" + wind->meteofile());
   noRecords = wind->noRecords_();
   ifstream* windparams_ptr = &windParamFile;
   windParamFile.close();
 }

 void concentration::prepare(double* posi, double relTime, double time, unsigned long id)
 {
   if ((nop > 0) && (relTime > currRelTime))
   {
    meanDrift[0] /= nop; meanDrift[1] /= nop; meanDrift[2] /= nop;
    for (int i = 0; i<nop; i++)
    {
	double sd = sqrt((position[i]-meanDrift[0])*(position[i]-meanDrift[0])+(position[ppr+i]-meanDrift[1])*(position[ppr+i]-meanDrift[1])+(position[2*ppr+i]-meanDrift[2])*(position[2*ppr+i]-meanDrift[2]));
	varPos += sd*sd;
	list<double>::iterator it = distrSigma.begin();
	
	while ( (it != distrSigma.end()) && (sd > *it) )
	{
	  it++;
	}
	distrSigma.insert(it,sd);
    }

    varPos=sqrt(varPos/nop);
    double h;
    if (nop > 100)
    {
      int m25=ceil(0.25*nop), m75=ceil(0.75*nop);
      list<double>::iterator it = distrSigma.begin();

      for (int i = 0; i < m25; i++)
	it++;

      double IQR = *it;

      for (int i = m25; i < m75; i++)
	it++;

      IQR = (*it - IQR)/1.34;

      if (varPos < IQR) 
        h = varPos;
      else
        h = IQR;
    }
    else
	h = varPos;

    h *= pow(nop,-1.0/7.0)*weight("bw");
    bandWidth.push_back(h);

    currRelTime = relTime;
    varPos = 0;
    nop = 0;
    distrSigma.clear();
    meanDrift[0] = 0; meanDrift[1] = 0; meanDrift[2] = 0;
   }
   else if ( distrSigma.size() == 0 )
   {
     currRelTime = relTime;
   }

   if ((abs(posi[0]) < xMax) && (abs(posi[1]) < yMax) && (abs(posi[2]) < zMax))
   {
    id = id - (id/ppr)*ppr;
    meanDrift[0] += posi[0]; meanDrift[1] += posi[1]; meanDrift[2] += posi[2];
    position[id] = posi[0]; position[ppr+id] = posi[1]; position[2*ppr+id] = posi[2];
    nop++;

    ofstream pfile(pfileName,ios::app);
    pfile << bandWidth.size() << "\t" << posi[0] << "\t" << posi[1] << "\t" << posi[2] << "\n" << endl;
    pfile.close();
   }
 }

 double concentration::weight(string type)
 {
   double w;

   if (type == "bw")
   {
     w = 720.0/675675; //beta value
     w *= 3*12006225/(4*acos(-1)*36); //inverse alpha squared value
     w *= 32*pow(acos(-1),1.5)/(3*(2+3));
     w = 0.85*pow(3*w,1.0/7.0);
   }
   else if (type == "kernel")
   {
     w = 0.25*3*315/(48*acos(-1));
   }

   return w;
 }

 int concentration::fac(int n)
 {
   int nfac=1;

   if (n==0)
	return nfac;

   for (int i = 1; i < n+1; i++)
	nfac *= i;

   return nfac;
 }

 double concentration::integralKernel(double zc, double h)
 {
   int a = 3; double intK=0, zLow, zUp;

   if (zc < h)
	zLow = -zc/h;
   else
	zLow = -1;

   if (zc+h > zMax)
	zUp = (zMax-zc)/h;
   else
	zUp = 1;

   for (int k = 0; k < a+2; k++)
     intK += pow(-1,k)*(pow(zUp,2*k+1)-pow(zLow,2*k+1))/((2*k+1)*fac(k)*fac(a+1-k));

   intK *= fac(a)*acos(-1);

   return intK;
 }

 void concentration::calculate(double Qp, double xMax, double zMax)
 {
    int Nx = 1, Ny = 1, Nz = 1; double dx=0, dy = 0, dz = 0;

    if ( Xmax > Xmin && resX > 0 )
    {
	Nx = ceil((Xmax-Xmin)/resX);
	if (Nx == 1)
	{
	  dx = Xmax-Xmin; Nx = 2;
	}
	else
	  dx = (Xmax-Xmin)/(Nx-1);
    }

    if ( Ymax > Ymin && resY > 0 )
    {
	Ny = ceil((Ymax-Ymin)/resY);
	if (Ny == 1)
	{
	  dy = Ymax-Ymin; Ny = 2;
	}
	else
	  dy = (Ymax-Ymin)/(Ny-1);
    }

    if ( Zmax > Zmin && resZ > 0 )
    {
	Nz = ceil((Zmax-Zmin)/resZ);
	if (Nz == 1)
	{
	  dz =  Zmax-Zmin; Nz = 2;
	}
	else
	  dz = (Zmax-Zmin)/(Nz-1);
    }

    double* c = new double[Nx*Ny*Nz](); 
    double x, y, z, Dzj, Dxj, dTd, Cda; 
    int ixLow, ixUp, iyLow, iyUp, izLow, izUp, bwIndex;
    ifstream pfile(pfileName); //,ios::app);
    while(!pfile.eof())
    {
	pfile >> bwIndex >> x >> y >> z;
	double h = bandWidth[bwIndex];
	if ( z < (Zmax+h) && z > (Zmin-h) && x < (Xmax+h) && x > (Xmin-h) && y < (Ymax+h) && y > (Ymin-h))
	{
	  ixLow = 0, iyLow = 0, izLow = 0, ixUp = 0, iyUp = 0, izUp = 0;
	  if (dx > 0) {ixLow = (x-h-Xmin)/dx; ixUp = ceil((x+h-Xmin)/dx);} if (ixLow < 0){ixLow=0;} if (ixUp > (Nx-1)){ixUp=Nx-1;}
	  if (dy > 0) {iyLow = (y-h-Ymin)/dy; iyUp = ceil((y+h-Ymin)/dy);} if (iyLow < 0){iyLow=0;} if (iyUp > (Ny-1)){iyUp=Ny-1;}
	  if (dz > 0) {izLow = (z-h-Zmin)/dz; izUp = ceil((z+h-Zmin)/dz);} if (izLow < 0){izLow=0;} if (izUp > (Nz-1)){izUp=Nz-1;}

	  for (int k = izLow; k < (izUp+1); k++)
	  {
	    Cda = integralKernel(Zmin+dz*k,h);
	    Dzj=(Zmin+dz*k-z)/h;
	    for (int j = ixLow; j < (ixUp+1); j++) // int i = ixLow; i < ixUp; i++)
	    {
	      Dxj=(Xmin+dx*j-x)/h;
	      for (int i = iyLow; i < (iyUp+1); i++) // int j = izLow; j < izUp; j++)
	      {
		dTd=(Ymin+dy*i-y)/h; dTd = Dxj*Dxj+Dzj*Dzj+dTd*dTd;
	  	if ( dTd < 1 )
	  	{
	    	  c[i+j*Ny+k*Nx*Ny] += pow(1-dTd,3)/(h*h*h*Cda);
	  	}
	      }
	    }
          }
	}
    }
    pfile.close();

    ofstream cfile(ofileName); double ci;
    for (unsigned i = 0; i < Nx*Ny*Nz; i++)
    {
	ci=c[i]*Qp;

	cfile << ci << setw(15) << endl;
    }
    delete[] c;
 }

 void concentration::reconstruct(int noProcs, string meteofile, int Nx, int Ny, int Nz, double dx, double dy, double dz, double xmin, double ymin, double zmin)
 {
   stringstream ss; string id; unsigned N = Nx*Ny*Nz; double* reconstrConc = new double[N]();
   for ( int j = 0; j < noProcs; j++)
   {
     ss << j; ss >> id;
     ifstream ifile("processor" + id + "/output/CONC_" + meteofile);
     for ( unsigned r=0; r < N; r++)
     {
        double ci; ifile >> ci;
        reconstrConc[r] += ci;
     }
   }

   ofstream ofile("output/CONC_" + meteofile); ofile.precision(6);
   int ixy, ix, iy, iz;
   for ( unsigned i = 0; i < N; i++)
   {
	ixy = (i % (Nx*Ny)); iy = (ixy % Ny); ix = ixy/Ny; iz = i/(Nx*Ny);
        ofile << xmin+ix*dx << setw(15) << ymin+iy*dy << setw(15) << zmin+iz*dz << setw(15) << reconstrConc[i];

        ofile << endl;
   }
   ofile.close();

   delete[] reconstrConc;
 }

