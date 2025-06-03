#!/bin/bash
#$ -q long
#$ -pe smp 24
#$ -N BFD


module use -a ~/privatemodules

module load beast267

# three different models, uncomment the one you want to run

beast -beagle -threads 24 inputsbfd_EWN.xml
#beast -beagle -threads 24 inputsbfd_N.xml
#beast -beagle -threads 24 inputsbfd_lineata.xml

echo "end"
