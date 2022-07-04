 #include <iostream>
 #include <fstream>
 #include <sstream>
 #include <cmath>
 #include "dictionary.h"
 #include "radiation.h"

int main()
{
   int noRecords=1;

   dictionary dict({"timeStep","measureFreq","noPartPerRel","sourceStrength","noProcs","meteofile"});
   int noprocs = dict.getfl(4);
   int ppr = ceil(dict.getfl(2)/noprocs);
   ppr = ppr/2 + (ppr % 2); ppr *= 2*noprocs;
   int mfreq = ceil(dict.getfl(1)/dict.getfl(0));
   double dt = dict.getfl(1)/mfreq;
   double Qp = dict.getfl(3)*dt/ppr;

   radiation dose(noRecords);
   for (int i = 0; i < noprocs; i++)
   {
     stringstream ss; ss << i; string procnr; ss >> procnr;
     ifstream posfile("processor" + procnr + "/output/Positions_" + dict.get(5));
     if (posfile.is_open())
     {
       cout << "Reading output processor " << i << endl;
       double ign, x, y, z, pos[3];
       bool no_eof=true;
       while (no_eof)
       {
          posfile >> ign >> x >> y >> z;
          if (!posfile.eof())
	  {
	    pos[0]=x; pos[1]=y; pos[2]=z;
	    dose.update(pos,noRecords,Qp);
	  }
          else
            no_eof=false;
       }
     }
     else
     {
	cout << "File 'Positions_" << dict.get(5) << "' not found in directory './output'." << endl;
     }
     posfile.close();
   }
   dose.write(dict.get(5));

   return 0;
}

