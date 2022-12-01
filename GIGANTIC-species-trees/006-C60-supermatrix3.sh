#!/bin/bash
# GIGANTIC species tree script

rm 5-input
for file in output/5-catsequences/*; do cat $file >> output/5-catsequences/5-input; printf '%s\n' >>  output/5-catsequences/5-input; done 

iqtree -s output/5-catsequences/5-input -m LG+C60 --prefix output/6-iqtree/speciestree --alrt 2000 -B 5000 -T 30 --threads-max 30
