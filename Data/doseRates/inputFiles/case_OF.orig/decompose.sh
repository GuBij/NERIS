#!/bin/bash -l 

exec=$VSC_DATA/C++/particleModel/simDose
meteoFile=meteo<day>.txt
option=Langevin

cd $VSC_SCRATCH/<folder>
k=$PBS_VNODENUM
mkdir processor$k
cp -r output processor$k/.
cp inputDict processor$k/.
cp $meteoFile processor$k/.
cd processor$k
$exec $k $option

