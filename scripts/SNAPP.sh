#!/bin/bash
#$ -q long
#$ -pe smp 48
#$ -N SNAPP


module use -a ~/privatemodules

module load beast277


echo "start" 

beast -beagle -threads 48 inputs/snapp.xml

echo "end"
