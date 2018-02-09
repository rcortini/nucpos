#!/bin/bash

for file in $(find . -name "*.dat"); do
  dist=$(echo $file | cut -d'_' -f1 | sed -e s,'./',, | tr '[:upper:]' '[:lower:]')
  mechanical_model=$(echo $file | cut -d'_' -f2)
  temperature=$(echo $file | cut -d'_' -f3 | sed -e s,'.dat',,)
  newfilename="$mechanical_model"_"$temperature".$dist
  cp $file $newfilename
done
