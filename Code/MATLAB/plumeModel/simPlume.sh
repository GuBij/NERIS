#!/bin/bash -l 

fileList=meteo/dateList;

cd $VSC_DATA/MATLAB/Scripts/plumeModel

k=$PBS_VNODENUM;
j=0;
while [ $j -le $k ] && IFS='' read -r line && [ -n "$line" ]; do
 j=`expr $j + 1`
done  < "$fileList"

./simPlume $line

# outputDir=output/$line;
# mkdir $outputDir

# Q=1;
# h=10;
# scheme=BM;
# nop=25;
# exec=../fluenceTools/rotConc;
# outputDir=$VSC_DATA/MATLAB/Output/NERIS/article/src60m_gauss/conc;

# for i in `seq 1 $nop`
# do
#  $exec meteo/meteo$line $i monitoringPoints $line
#  pointFieldConc $line $Q $h $scheme meteo
#  mv $outputDir/pointFieldConc.out.txt $outputDir/conc"$line"_p$i
# done

