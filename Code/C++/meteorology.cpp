#include "Uforest.h"
#include "UopenField.h"
#include "KTaylorOF.h"
#include "Kforest.h"
#include "KopenField.h"
#include "KTaylorForest.h"

 meteorology::meteorology ()
 {}

 meteorology* meteorology::New(bool LangevinModel)
 {
   dictionary dict({"meteofile","terrainModel","stackHeight","diameter","flowRate","gasTemp"});

   if ( dict.get(1) == "openField" )
	return new openField(dict, LangevinModel);
   else if ( dict.get(1) == "forest" )
        return new Uforest(dict, LangevinModel);
   else
	return nullptr;
 }

 meteorology* meteorology::New(string type, bool LangevinModel)
 {
   dictionary dict({"meteofile","terrainModel","stackHeight","timeStep","diameter","flowRate","gasTemp"});

   string modelName = dict.get(1);
   if (type == "windModel")
   {
	if (modelName == "openField")
	   return new UopenField(dict);
	else if (modelName == "forest")
	   return new Uforest(dict);
	else
	   return nullptr;
   }
   else if (type == "eddyDiffModel")
   {
        if (modelName == "openField" && LangevinModel)
	   return new KTaylorOF(dict);
	else if (modelName == "openField")
           return new KopenField(dict);
	else if (modelName == "forest" && LangevinModel)
	   return new KTaylorForest(dict);
        else if (modelName == "forest")
           return new Kforest(dict);
	else
	   return nullptr;
   }
   else
	error msg("inputDict",type);
 }

