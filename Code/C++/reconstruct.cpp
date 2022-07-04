 #include <iostream>
 #include <cmath>
 #include <string>
 #include "radiation.h"
 #include "concentration.h"
 #include "dictionary.h"

 int main(int argc, char **argv)
 {
   string option(argv[argc-1]);
   dictionary dict({"meteofile","noProcs","xMax","zMax"});

   if (option == "dose")
	radiation::reconstruct(int(dict.getfl(1)),dict.get(0));
   else if (option == "conc")
   {
	dictionary sampleDict(dict,"sampleConc",{"Xmin","Xmax","Ymin","Ymax","Zmin","Zmax","resX","resY","resZ"});

	double Xmin = sampleDict.getfl(0), Xmax = sampleDict.getfl(1), Ymin = sampleDict.getfl(2), Ymax = sampleDict.getfl(3), Zmin = sampleDict.getfl(4), Zmax = sampleDict.getfl(5);
	double resX = sampleDict.getfl(6), resY = sampleDict.getfl(7), resZ = sampleDict.getfl(8);
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

	concentration::reconstruct(int(dict.getfl(1)),dict.get(0),Nx,Ny,Nz,dx,dy,dz,Xmin,Ymin,Zmin);
   }
   else
   {
	cout << "\tOption " << argv[argc-1] << " not implemented. Valid options are \n"
	     << "\t  dose\n"
	     << "\t  conc\n" << endl;
   }
   return 0;
 }

