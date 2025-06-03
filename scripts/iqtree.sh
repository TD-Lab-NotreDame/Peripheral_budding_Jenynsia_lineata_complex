#!/bin/bash
#$ -q long
#$ -pe smp 4
#$ -N iqtree

conda activate environment_name  # Replace with your actual conda environment name

#First run IQtree to generate a file that will pass filters for only containing variant sites when using an ASC adjustment \
#No outgroup\
iqtree -s inputs/filtered_final_83_LD_pruned.recode.min4.phy -st DNA -m GTR+G4+F+ASC -bb 1000 -alrt 1000
#With an outgroup \
iqtree -s inputs/filtered_final_94_LD_pruned.recode.min4.phy -st DNA -m GTR+G4+F+ASC -bb 1000 -alrt 1000
#Use the output files that end in ".varsites.phy" for further analysis \
#No outgroup iqtree\
wait
iqtree -s filtered_final_83_LD_pruned.recode.min4.phy.varsites.phy -st DNA -m GTR+G4+F+ASC -bb 1000 -alrt 1000

echo "end"
