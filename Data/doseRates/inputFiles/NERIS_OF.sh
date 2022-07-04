#!/bin/bash

  while IFS='' read -r line; do
  option=Langevin
        if [[ ( -n $line )  && ( $line != *%* ) ]];then
          dirName=${line//[ ,   ]/}
          cp -r case_OF.orig $dirName
          cd $dirName
          echo preparing directory $dirName
          cp meteo/meteo$dirName* .
          sed -i "s/<day>/$dirName/" inputDict
          sed -i "s/<nameJob>/$dirName/" script.pbs
          sed -i "s/<folder>/$dirName/" script.pbs
		  sed -i "s/<folder>/$dirName/" decompose.sh
		  sed -i "s/<day>/$dirName/" decompose.sh
          $VSC_DATA/C++/particleModel/genMeteo $option
          cd ../.
        fi
  done < "$1"

