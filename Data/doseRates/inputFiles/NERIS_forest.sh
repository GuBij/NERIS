#!/bin/bash -l

if [[ ( $# -eq 2 ) && ( $2 == -qsub  ) ]]; then
  while IFS='' read -r line; do
        if [[ ( -n $line )  && ( $line != *%* ) ]];then
	  allvalues=()
#          dirName=${line//[ ,   ]/}
	  allvalues+=( $line )
	  dirName=${allvalues[0]}
          cd $dirName
          qsub script.pbs
          cd ..
        fi
  done < "$1"
elif [[ ( $# -eq 2 ) && ( $2 == -recon  ) ]]; then
  while IFS='' read -r line; do
        if [[ ( -n $line )  && ( $line != *%* ) ]];then
	  allvalues=()
#          dirName=${line//[ ,   ]/}
          allvalues+=( $line )
          dirName=${allvalues[0]}
          cd $dirName
          $VSC_DATA/C++/particleModel/reconstruct dose
	  cd output
	  mv Dose_meteo$dirName".txt" $dirName"17.txt"
	  echo Dose $dirName reconstructed
          cd ../..
        fi
  done < "$1"
elif [[ ( $# -eq 2 ) && ( $2 == -cp  ) ]]; then
  while IFS='' read -r line; do
        if [[ ( -n $line )  && ( $line != *%* ) ]];then
	  allvalues=()
#          dirName=${line//[ ,   ]/}
	  allvalues+=( $line )
	  dirName=${allvalues[0]}
          cd $dirName/output
          cp $dirName"17.txt" $VSC_DATA/MATLAB/Output/NERIS/article/PM/src60m/forest/.
	  echo Copied $dirName output
          cd ../..
        fi
  done < "$1"
elif [[ $# -eq 1 ]]; then
  option=Langevin
  while IFS='' read -r line; do
        if [[ ( -n $line )  && ( $line != *%* ) ]];then
	  allvalues=()
#          dirName=${line//[ ,   ]/}
	  allvalues+=( $line )
	  dirName=${allvalues[0]}
	  LAI=${allvalues[1]}
          cp -r case_forest.orig $dirName
          cd $dirName
	  echo ============================
          echo preparing directory $dirName
	  echo ============================
          cp meteo/meteo$dirName* .
          sed -i "s/<day>/$dirName/" inputDict
		  sed -i "s/<LAI>/$LAI/" inputDict
          sed -i "s/<nameJob>/$dirName/" script.pbs
          sed -i "s/<folder>/$dirName/" script.pbs
		  sed -i "s/<folder>/$dirName/" decompose.sh
		  sed -i "s/<day>/$dirName/" decompose.sh
          $VSC_SCRATCH/PM/genMeteo $option
          cd ../. 
        fi
  done < "$1"
fi
