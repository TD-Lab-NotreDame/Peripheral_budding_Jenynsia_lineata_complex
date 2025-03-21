##Code for genomic analyses carreid out in terminal

### analyses using the 83 individuals dataset

#################### ddocent to trim, align, map, and call snps ########################################################
##################################################################################################################
# denovo approach
# This was run in a conda environemnt 
# ddocent needs to be installed  https://ddocent.com/bioconda/
# file ddocent_config_denovo.sh needs to be in working directory
# fastq files must be in the working directory (excluding J. onscura samples)

dDocent ddocent_config_denovo.sh

############################ filtering vcftools ##################################################################
##################################################################################################################
# vcftools needs to be installed in conda

#calculate missingess for each individual on raw ddocent output
vcftools --vcf TotalRawSNPs.vcf --missing-indv
#from output make a list of what samples have more than 50% missing data
mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
#removed from ddocent filtered vcf file
vcftools --vcf ddocent_Final.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out filter1

##look at population levels of missing loci
# popmap83 must be in working directory
mawk '$2 == "lineataE"' popmap83 > 1.keep && mawk '$2 == "lineataW"' popmap83 > 2.keep && mawk '$2 == "AFlineata"' popmap83 > 3.keep \
&& mawk '$2 == "onca"' popmap83 > 4.keep && mawk '$2 == "luxata"' popmap83 > 5.keep && mawk '$2 == "darwinii"' popmap83 > 6.keep
# filter missing loci by population
vcftools --vcf filter1.recode.vcf --keep 1.keep --missing-site --out 1
vcftools --vcf filter1.recode.vcf --keep 2.keep --missing-site --out 2 
vcftools --vcf filter1.recode.vcf --keep 3.keep --missing-site --out 3
vcftools --vcf filter1.recode.vcf --keep 4.keep --missing-site --out 4
vcftools --vcf filter1.recode.vcf --keep 5.keep --missing-site --out 5
vcftools --vcf filter1.recode.vcf --keep 6.keep --missing-site --out 6 

#remove loci with more than 10% missing data
cat 1.lmiss 2.lmiss 3.lmiss 4.lmiss 5.lmiss 6.lmiss 7 | mawk '!/CHR/' | mawk '$6 > 0.1' | cut -f1,2 >> badloci
vcftools --vcf filter1.recode.vcf --exclude-positions badloci --recode --recode-INFO-all --out filter2

#appluy RAD sequencing filters
#vcffilter must be installed in conda
vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" filter2.recode.vcf > filter3.vcf
vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s filter3.vcf > Jfilter4.vcf
vcffilter -f "MQM / MQMR > 0.9 & MQM / MQMR < 1.05" filter4.vcf > filter5.vcf
##important for de novo
vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s filter5.vcf > filter6.vcf

#final filters
VCF_IN=filter6.vcf
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

#LD pruning for 1 snp per 500 bp
vcftools --vcf filtered_final_83.recode.vcf \
 --thin 500 --recode --out filtered_final_83_LD_pruned



########################### file conversion ######################################################################
##################################################################################################################
# use vcf2phylip (download to conda) to convert to nexus, fasta, and phylip formats

python software/vcf2phylip.py -i filtered_final_83_LD_pruned.recode.vcf -f -n -b

########################### PCA ######################################################################
##################################################################################################################
#ipyrad must be installed in conda
#compress vcf file
bgzip filtered_final_83.recode.vcf > filtered_final_83.recode.vcf.gz
# tabix index the compressed VCF (creates .vcf.gz.tbi)
tabix filtered_final_83.recode.vcf.gz

# run python file ipyrad_converter_pca
# run python file ipyrad_pca.py
##outputs PCA figure

########################### admixture ######################################################################
##################################################################################################################
#PLINK needs to be installed in conda

#generate PLINK format files
VCF=filtered_final_83_LD_pruned.recode.vcf

plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# --vcf-half-call missing --make-bed \
#--indep-pairwise 50 10 0.2 \
--out plink

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

#outputs CV errors and files for plotting each K in r

#generate admixture plots in r

################################# RAXML ###############################################################
##################################################################################################################
#RAXML must be downloaded to conda

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 1000 -s filtered_final_83_LD_pruned.recode.min4.phy -n T20

########################### SVDquartets ######################################################################
##################################################################################################################
# PAUP needs to be installed and opened
# svd.nexus needs to be in the working directory
# taxpartitions.txt needs to be in the working directory
#open file
exe svd.nexus
#run
svdq taxpartition=fish showScores=no seed=1234568 bootstrap nreps=1000 treeFile=svd.tre;
#save
savetree

### analyses using the 30 individuals dataset
# need snapp_remove.txt in working directory
vcftools --vcf filter6.vcf --remove snapp_remove.txt --recode --recode-INFO-all --out filtered_final_30

VCF_IN=filtered_final_30.recode.vcf
VCF_OUT=filtered_final_30_nomiss_LD_pruned

MAF=0.1
MISS=1
QUAL=30
MIN_DEPTH=12
MAX_DEPTH=63
MIN_A_COUNT=3 

vcftools --vcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-alleles 2 --max-alleles 2 \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --mac $MIN_A_COUNT --thin 500 \
--recode-INFO-all --recode --out $VCF_OUT

mawk '!/#/' Jenynsia_rad_denovo_all_filtered_nomiss_LD_pruned.recode.vcf | wc -l

# Randomly select 1000 lines from a VCF file and keep them in original order
grep '^#' filtered_final_30_nomiss_LD_pruned.recode.vcf > filtered_final_30_nomiss_LD_pruned_sub1000.vcf  # Copy header lines to output

grep -v '^#' filtered_final_30_nomiss_LD_pruned.recode.vcf | shuf | head -n 1000 | sort -n | tail -n 1000 >> filtered_final_30_nomiss_LD_pruned_sub1000.vcf

# file conversion
python software/vcf2phylip.py -i filtered_final_30_nomiss_LD_pruned_sub1000.vcf-f -n -b
#binary nexus used for SNAPP analysis

################################# SNAPP and BFD ###############################################################
##################################################################################################################
#Beast 2 needs to be downloaded
# to run in beast
# xml files need to be in working directory
beast -threads 4 file.xml

################################# DELINEATE ###############################################################
##################################################################################################################
#delineate needs to be downloaded to conda environment
#delin_guide_tree_final.nex and delin_input_final.tsv must be in working directory
#guide tree is consensus tree from the 2 SNAPP runs, generated in tree annotater 

delineate-estimate partitions --tree-file delin_guide_tree_final.nex --constraints delin_input_final.tsv




### analyses using the 94 individuals dataset
#follow same filtering steps as above, but include J. obscura fastq files in working directory before running ddocent.

dDocent ddocent_config_denovo.sh

############################ filtering vcftools ##################################################################
##################################################################################################################
# vcftools needs to be installed in conda

#calculate missingess for each individual on raw ddocent output
vcftools --vcf TotalRawSNPs.vcf --missing-indv
#from output make a list of what samples have more than 50% missing data
mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
#removed from ddocent filtered vcf file
vcftools --vcf ddocent_Final.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out filter1

##look at population levels of missing loci
# popmap83 must be in working directory
mawk '$2 == "lineataE"' popmap94 > 1.keep && mawk '$2 == "lineataW"' popmap94 > 2.keep && mawk '$2 == "AFlineata"' popmap94 > 3.keep \
&& mawk '$2 == "onca"' popmap94 > 4.keep && mawk '$2 == "luxata"' popmap94 > 5.keep && mawk '$2 == "darwinii"' popmap94 > 6.keep && mawk '$2 == "obscura"' popmap94 > 7.keep
# filter missing loci by population
vcftools --vcf filter1.recode.vcf --keep 1.keep --missing-site --out 1
vcftools --vcf filter1.recode.vcf --keep 2.keep --missing-site --out 2 
vcftools --vcf filter1.recode.vcf --keep 3.keep --missing-site --out 3
vcftools --vcf filter1.recode.vcf --keep 4.keep --missing-site --out 4
vcftools --vcf filter1.recode.vcf --keep 5.keep --missing-site --out 5
vcftools --vcf filter1.recode.vcf --keep 6.keep --missing-site --out 6 
vcftools --vcf filter1.recode.vcf --keep 6.keep --missing-site --out 7 

#remove loci with more than 10% missing data
cat 1.lmiss 2.lmiss 3.lmiss 4.lmiss 5.lmiss 6.lmiss 7 | mawk '!/CHR/' | mawk '$6 > 0.1' | cut -f1,2 >> badloci
vcftools --vcf filter1.recode.vcf --exclude-positions badloci --recode --recode-INFO-all --out filter2

#appluy RAD sequencing filters
#vcffilter must be installed in conda
vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" filter2.recode.vcf > filter3.vcf
vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s filter3.vcf > Jfilter4.vcf
vcffilter -f "MQM / MQMR > 0.9 & MQM / MQMR < 1.05" filter4.vcf > filter5.vcf
##important for de novo
vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s filter5.vcf > filter6.vcf

#final filters
VCF_IN=filter6.vcf
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

#LD pruning for 1 snp per 500 bp
vcftools --vcf filtered_final_94.recode.vcf \
 --thin 500 --recode --out filtered_final_94_LD_pruned

########################### file conversion ######################################################################
##################################################################################################################
# use vcf2phylip (download to conda) to convert to nexus, fasta, and phylip formats

python software/vcf2phylip.py -i filtered_final_94_LD_pruned.recode.vcf -f -n -b

########################### PCA ######################################################################
##################################################################################################################
#ipyrad must be installed in conda
#compress vcf file
bgzip filtered_final_94.recode.vcf > filtered_final_94.recode.vcf.gz
# tabix index the compressed VCF (creates .vcf.gz.tbi)
tabix filtered_final_94.recode.vcf.gz
################################# Treemix ###############################################################
##################################################################################################################

#run python file ipyrad_converter_treemix.py
#run python file ipyrad_treemix.py

################################# DSUITE ###############################################################
##################################################################################################################
# download dsuite https://github.com/millanek/Dsuite
# svd.newick.txt must be in working directory, output from SVDquartets
##step 1
Dsuite Dtrios -o Jenynsia_rad_denovo_dsuite -t svd.newick.txt filtered_final_94_LD_pruned..recode.vcf.gz sets.txt 
##step 2
Dsuite Fbranch svd.newick.txt dsuite_tree.txt > output
#outputs dsuite.fb.txt

#to plot
#dtools.py is part of the dsuite software
python3 dtools.py dsuite.fb.txt svd.newick.txt

