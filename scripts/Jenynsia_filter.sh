#!/bin/bash
#$ -q long
#$ -pe smp 4
#$ -N filtering

conda activate environment_name  # Replace with your actual conda environment name

echo "start"
########for 83 individuals dataset ###########################################
#calculate missingess for each individual on raw ddocent output
vcftools --vcf TotalRawSNPs_83.vcf --missing-indv
wait
#from output make a list of what samples have more than 50% missing data
mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP83.indv
#removed from ddocent filtered vcf file
wait
vcftools --vcf ddocent_Final_83.recode.vcf --remove lowDP83.indv --recode --recode-INFO-all --out filter1_83
wait
##look at population levels of missing loci
# popmap83 must be in working directory
mawk '$2 == "lineataE"' popmap83 > 1.keep && mawk '$2 == "lineataW"' popmap83 > 2.keep && mawk '$2 == "AFlineata"' popmap83 > 3.keep \
&& mawk '$2 == "onca"' popmap83 > 4.keep && mawk '$2 == "luxata"' popmap83 > 5.keep && mawk '$2 == "darwinii"' popmap83 > 6.keep
wait
# filter missing loci by population
vcftools --vcf filter1_83.recode.vcf --keep 1.keep --missing-site --out 1
vcftools --vcf filter1_83.recode.vcf --keep 2.keep --missing-site --out 2 
vcftools --vcf filter1_83.recode.vcf --keep 3.keep --missing-site --out 3
vcftools --vcf filter1_83.recode.vcf --keep 4.keep --missing-site --out 4
vcftools --vcf filter1_83.recode.vcf --keep 5.keep --missing-site --out 5
vcftools --vcf filter1_83.recode.vcf --keep 6.keep --missing-site --out 6 
wait
#remove loci with more than 10% missing data
cat 1.lmiss 2.lmiss 3.lmiss 4.lmiss 5.lmiss 6.lmiss 7 | mawk '!/CHR/' | mawk '$6 > 0.1' | cut -f1,2 >> badloci
wait
vcftools --vcf filter1_83.recode.vcf --exclude-positions badloci --recode --recode-INFO-all --out filter2_83
wait
#appluy RAD sequencing filters
#vcffilter must be installed in conda
vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" filter2_83.recode.vcf > filter3_83.vcf
wait
vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s filter3_83.vcf > filter4_83.vcf
wait
vcffilter -f "MQM / MQMR > 0.9 & MQM / MQMR < 1.05" filter4_83.vcf > filter5_83.vcf
wait
##important for de novo
vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s filter5.vcf > filter6_83.vcf
wait
#final filters
VCF_IN=filter6_83.vcf
VCF_OUT=filtered_final_83

MISS=0.9
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=58
MIN_A_COUNT=3

vcftools --vcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-alleles 2 --max-alleles 2 \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --mac $MIN_A_COUNT \
--recode-INFO-all --recode --out $VCF_OUT
wait
#LD pruning for 1 snp per 500 bp
vcftools --vcf filtered_final_83.recode.vcf \
 --thin 500 --recode --out inputs/filtered_final_83_LD_pruned
wait
### analyses using the 30 individuals dataset
# need snapp_remove.txt in working directory
vcftools --vcf filter6_83.vcf --remove inputs/snapp_remove.txt --recode --recode-INFO-all --out filtered_final_30

VCF_IN=filtered_final_30.recode.vcf
VCF_OUT=filtered_final_30_nomiss_LD_pruned

MAF=0.1
MISS=1
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=58
MIN_A_COUNT=3 

vcftools --vcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-alleles 2 --max-alleles 2 \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --mac $MIN_A_COUNT --thin 500 \
--recode-INFO-all --recode --out $VCF_OUT

mawk '!/#/' filtered_final_30_nomiss_LD_pruned.recode.vcf | wc -l

# Randomly select 1000 lines from a VCF file and keep them in original order
grep '^#' filtered_final_30_nomiss_LD_pruned.recode.vcf > inputs/filtered_final_30_nomiss_LD_pruned_sub1000.vcf  # Copy header lines to output

grep -v '^#' filtered_final_30_nomiss_LD_pruned.recode.vcf | shuf | head -n 1000 | sort -n | tail -n 1000 >> inputs/filtered_final_30_nomiss_LD_pruned_sub1000.vcf






#########for 94 individuals dataset ###########################################
##uncomment to run for 94 individuals
#vcftools --vcf TotalRawSNPs_94.vcf --missing-indv
#wait
##from output make a list of what samples have more than 50% missing data
#mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP_94.indv
#wait
##removed from ddocent filtered vcf file
#vcftools --vcf ddocent_Final_94.recode.vcf --remove lowDP_94.indv --recode --recode-INFO-all --out filter1_94
#wait
###look at population levels of missing loci
## popmap83 must be in working directory
#mawk '$2 == "lineataE"' popmap94 > 1.keep && mawk '$2 == "lineataW"' popmap94 > 2.keep && mawk '$2 == "AFlineata"' popmap94 > 3.keep \
#&& mawk '$2 == "onca"' popmap94 > 4.keep && mawk '$2 == "luxata"' popmap94 > 5.keep && mawk '$2 == "darwinii"' popmap94 > 6.keep && mawk '$2 == "obscura"' popmap94 > 7.keep
#wait
## filter missing loci by population
#vcftools --vcf filter1_94.recode.vcf --keep 1.keep --missing-site --out 1
#vcftools --vcf filter1_94.recode.vcf --keep 2.keep --missing-site --out 2 
#vcftools --vcf filter1_94.recode.vcf --keep 3.keep --missing-site --out 3
#vcftools --vcf filter1_94.recode.vcf --keep 4.keep --missing-site --out 4
#vcftools --vcf filter1_94.recode.vcf --keep 5.keep --missing-site --out 5
#vcftools --vcf filter1_94.recode.vcf --keep 6.keep --missing-site --out 6 
#vcftools --vcf filter1_94.recode.vcf --keep 6.keep --missing-site --out 7 
#wait
##remove loci with more than 10% missing data
#cat 1.lmiss 2.lmiss 3.lmiss 4.lmiss 5.lmiss 6.lmiss 7 | mawk '!/CHR/' | mawk '$6 > 0.1' | cut -f1,2 >> badloci
#wait
#vcftools --vcf filter1_94.recode.vcf --exclude-positions badloci --recode --recode-INFO-all --out filter2_94
#wait
##appluy RAD sequencing filters
##vcffilter must be installed in conda
#vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" filter2_94.recode.vcf > filter3_94.vcf
#wait
#vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s filter3_94.vcf > filter4_94.vcf
#wait
#vcffilter -f "MQM / MQMR > 0.9 & MQM / MQMR < 1.05" filter4_94.vcf > filter5_94.vcf
#wait
###important for de novo
#vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s filter5_94.vcf > filter6_94.vcf
#wait
##final filters
#VCF_IN=filter6.vcf
#VCF_OUT=filtered_final_94
#
#MISS=0.9
#QUAL=30
#MIN_DEPTH=10
#MAX_DEPTH=58
#MIN_A_COUNT=3
#
#vcftools --vcf $VCF_IN \
#--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-alleles 2 --max-alleles 2 \
#--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --mac $MIN_A_COUNT \
#--recode-INFO-all --recode --out $VCF_OUT
#
##LD pruning for 1 snp per 500 bp
#vcftools --vcf filtered_final_94.recode.vcf \
# --thin 500 --recode --out inputs/filtered_final_94_LD_pruned
#
#







echo "end"
