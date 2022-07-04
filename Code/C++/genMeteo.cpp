#include <iostream>
#include <fstream>
#include <memory>
#include "meteorology.h"

 int main(int argc, char **argv)
 {

   bool LangevinModel = false;
   if ((argc > 1) && (string(argv[1]) == "Langevin"))
   {
      LangevinModel = true;
      cout << "\nLangevin model selected.\n" << endl;
   }
   else if (argc > 1)
   {
      cout << "Option " << argv[1] << " not implemented." << endl;
   }

   unique_ptr<meteorology> metModel(meteorology::New(LangevinModel));

   if (metModel == nullptr)
   {
	cout << "\nError: entry for 'terrainModel' in inputDict is possibly erroneous.\n" << endl; 
	error msg();
   }

   ifstream meteofile(metModel->meteofile());
   if (meteofile.is_open())
   {
    double meteoData[7]; cout << "noRecords: " << metModel->noRecords_() << endl;
    for(int i = 1; i < metModel->noRecords_(); i++)
     {
        meteofile >> meteoData[0] >> meteoData[1] >> meteoData[2] >> meteoData[3] >> meteoData[4] >> meteoData[5] >> meteoData[6];
        metModel->calcParams(i-1,meteoData);
     }
   }
   else 
   {
     error msg(metModel->meteofile());
   }

   return 0;
 }
