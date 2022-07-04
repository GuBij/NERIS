#!/bin/bash -l

 if [[ $1 == -genMeteo ]]; then
  $VSC_DATA/C++/particleModel/genMeteo
 elif [[ ( $# -eq 2 ) && ( $1 == -recon ) ]]; then
  $VSC_DATA/C++/particleModel/reconstruct dose
  dirName=$2
  cd output
  mv Dose_meteo$dirName".txt" $dirName"17.txt"
  echo Dose $dirName reconstructed
  cd ../..
 elif
  echo Option -recon requires extra input or
  echo option not implemented
 fi

