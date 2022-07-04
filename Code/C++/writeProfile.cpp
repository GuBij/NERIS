
 #include "Kforest.h"

 int main()
 {
   dictionary dict({"meteofile","terrainModel","stackHeight","timeStep"});

   Kforest Uf(dict);

   ifstream* windPtr = new ifstream("output/PARAMS_" + Uf.name()+ "_" + Uf.meteofile());
   Uf.updateParams(windPtr);

   Uf.writeProfile();
   delete windPtr;

   return 0;
 }

