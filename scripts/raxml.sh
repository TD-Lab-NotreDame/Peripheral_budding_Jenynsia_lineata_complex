#!/bin/bash
#$ -q long
#$ -pe smp 4
#$ -N raxml

conda activate environment_name  # Replace with your actual conda environment name

#No outgroup RAXML \
raxmlHPC -f a -m ASC_GTRGAMMA --asc-corr=lewis -p 12345 -x 12345 -# 1000 -s filtered_final_83_LD_pruned.recode.min4.phy.varsites.phy -n T1
#With an outgroup RAXML \
raxmlHPC -f a -m ASC_GTRGAMMA --asc-corr=lewis -p 12345 -x 12345 -# 1000 -s filtered_final_94_LD_pruned.recode.min4.phy.varsites.phy -n T2

echo "end"
