#!/bin/bash
# GIGANTIC species tree script

mkdir output/1-list output/2-split output/3-mafft output/4-clipkit output/5-catsequences output/6-iqtree

ls output/eachbusco_X_allspecies/ | awk -F "-" '{print $NF}'  > output/1-list/busco-list
grep -h '>' output/eachbusco_X_allspecies/* | cut -d'-' -f6,7 | sort | uniq > output/1-list/species-list

var=$( ls output/eachbusco_X_allspecies/ | awk 'BEGIN{FS=OFS="-"}{NF--; print}'| sort | uniq)

split -n l/10 output/1-list/busco-list output/2-split/chunk-

for file in output/2-split/chunk-*
do
  while read p
  do
    linsi --quiet output/eachbusco_X_allspecies/"$var"-"$p" > output/3-mafft/mafft-$p #&&
    clipkit output/3-mafft/mafft-$p -o output/4-clipkit/clipkit-$p
  done < $file&
done
