#!/bin/bash
#$ -q long
#$ -pe smp 4
#$ -N dsuite

conda activate environment_name  # Replace with your actual conda environment name

Dsuite Dtrios -o Jenynsia_rad_denovo_dsuite -t inputs/svd.newick.txt inputs/filtered_final_94_LD_pruned.recode.vcf.gz inputs/sets.txt 
##step 2
Dsuite Fbranch svd.newick.txt dsuite_tree.txt > output
#outputs dsuite.fb.txt




echo "end"
