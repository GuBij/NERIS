#!/bin/bash -l 

exec=$VSC_DATA/C++/particleModel/simDose;
meteofile=meteo<day>.txt

cd $VSC_SCRATCH/NERIS_2/<folder>
k=$PBS_VNODENUM
mkdir processor$k
cp -r output processor$k/.
cp inputDict processor$k/.
cp $meteofile processor$k/.
cd processor$k
$exec $k

