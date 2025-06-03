#!/bin/bash
#$ -q long
#$ -pe smp 4
#$ -N delineate

conda activate environment_name  # Replace with your actual conda environment name

delineate-estimate partitions --tree-file inputs/delin_guide_tree_final.nex --constraints inputs/delin_input_final.tsv

echo "end"
