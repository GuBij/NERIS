#include <iostream>
#include <iomanip>
#include <sstream>
#include "particle.h"
#include "radiation.h"

 int main(int argc, char **argv)
 {
   int proc = 0;
   bool LangevinModel = false;
   if (argc > 1)
   {
     for ( int i = 1; i < argc; i++)
     {
        if ( string(argv[i]) == "Langevin")
	{
	  LangevinModel = true;
	  cout << "\nLangevin model selected.\n" << endl;
	}
	else
	{
          stringstream ss(argv[i]); ss >> proc;
	}
     }
   }

   meteorology* windModel(meteorology::New("windModel", LangevinModel));
   meteorology* eddyDiffModel(meteorology::New("eddyDiffModel", LangevinModel));

   dictionary dict({"stackHeight","timeStep","measureFreq","noPartPerRel","xMax","yMax","zMax","sourceStrength","noProcs"});
   double dictParamValues[9];
   for ( int i = 0; i < 9; i++ )
	dictParamValues[i]=dict.getfl(i);

   cout << "\n domain: [-" << dictParamValues[4] << " m, " << dictParamValues[4] << " m ] x [-" << dictParamValues[5] << " m, " << dictParamValues[5] << " m ] x [0 m, " << dictParamValues[6] << " m ]" << endl;
   cout << "\n source strength: " << dictParamValues[7] << " Bq/h " << endl; 
   cout << "\n stack height: " << dictParamValues[0] << " m" << endl;
   cout << "\n " << dictParamValues[3] << " particles per time step are released" << endl;

   dictParamValues[3] = ceil(dictParamValues[3]/dictParamValues[8]);

   radiation dose(windModel->noRecords_()-1);

   int i=0;
   int warmUpPeriod=0; 
   bool proceed(true);
   while ( proceed )
   {
      particle pi(proc,i,dictParamValues,windModel,eddyDiffModel,warmUpPeriod,LangevinModel);
      while ( pi.inside() && pi.alive() )
      {
	  pi.advance(); 
	  double loc[3]; 

	  if ( pi.measurement() )
	  {
	     loc[0]=pi.coord()[0]; loc[1]=pi.coord()[1]; loc[2]=pi.coord()[2]; dose.update(loc,pi.time()/dictParamValues[2],pi.particleCharge()); 
             loc[0]=pi.coord()[3]; loc[1]=pi.coord()[4]; loc[2]=pi.coord()[5]; dose.update(loc,pi.time()/dictParamValues[2],pi.particleCharge());
	  }
      }
      pi.close();
      proceed = !pi.eos(); 
      i++;
   }
   cout << "\n warming up period: " << warmUpPeriod*dictParamValues[1] << " s \n" << endl;

   dose.write(windModel->meteofile());

   delete windModel; delete eddyDiffModel;

   return 0;
 }
