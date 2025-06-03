#!/bin/bash
#$ -q long
#$ -pe smp 4
#$ -N ddocent

conda activate environment_name  # Replace with your actual conda environment name

dDocent scripts/ddocent_config_denovo.sh

echo "end"
