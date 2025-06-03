#!/bin/bash
#$ -q long
#$ -pe smp 4
#$ -N admixture

conda activate environment_name  # Replace with your actual conda environment name

#generate PLINK format files
VCF=inputs/filtered_final_83_LD_pruned.recode.vcf

plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# --vcf-half-call missing --make-bed \
--out plink
wait
# start setting up for for admixture
# https://speciationgenomics.github.io/ADMIXTURE/
awk '{$1="0";print $0}' plink.bim > plink.bim.tmp
mv plink.bim.tmp plink.bim
wait
FILE=plink

#fun K 2-10
for i in {2..10}
do
 admixture --cv $FILE.bed $i > log${i}.out
done

echo "end"
