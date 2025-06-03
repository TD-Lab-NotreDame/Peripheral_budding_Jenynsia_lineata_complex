The files here are associated with "Peripheral budding following range expansion explains diversity and distribution of one-sided livebearing fish". Preprint available at https://doi.org/10.22541/au.174255499.93236941/v1

The scripts folder contains the scripts used in analyses, the inputs folder contains inputs for genomic analyses, the data folder contains data files used for *.R files.

The code in the *.sh files was written for use on an HPC sungrid engine system.

Phylip and nexus files are also included from phylogenetic analyses.

Genomic analyses:

**dDocent**

This code was run in a conda environemnt \
dDocent needs to be installed  https://ddocent.com/bioconda/ \
fastq files must be in the working directory \
For the 83 and 30 individuals data set do not have J. obscura samples in the working directory \
For the 94 individuals dataset, include all fatsq files \
Run \
```scripts/ddocent.sh``` \
Rename output to TotalRawSNPs_83.vcf and ddocent_Final_83.recode.vcf after running on 83 individuals \
Rename output to TotalRawSNPs_94.vcf and ddocent_Final_94.recode.vcf after running on 94 individuals 

**Filtering VCF files**

vcftools and vcffilter need to be installed in conda \
run ```scripts/Jenynsia_filter.sh```


**File conversion**

download vcf2phylip https://github.com/edgardomortiz/vcf2phylip to convert to nexus, fasta, and phylip formats\
run \
for 83 individuals dataset \
```python vcf2phylip.py -i inputsfiltered_final_83_LD_pruned.recode.vcf -f -n -b```\
for 30 individuals dataset \
```python vcf2phylip.py -i inputs/filtered_final_30_nomiss_LD_pruned_sub1000.vcf -f -n -b```\
for 94 individuals dataset \
```python vcf2phylip.py -i inputs/filtered_final_94_LD_pruned.recode.vcf -f -n -b```

**PCA**

ipyrad needs to be installed in conda \
compress vcf file \
```bgzip inputs/filtered_final_83.recode.vcf > inputs/filtered_final_83.recode.vcf.gz``` \
index the compressed file \
```tabix inputs/filtered_final_83.recode.vcf.gz``` \
then run the python file ```ipyrad_converter_pca``` \
followed by ```ipyrad_pca.py```\
A PCA figure will be output

**Admixture**

PLINK needs to be installed in conda \
run ```scripts/Jenynsia_admixture.sh``` \
Outputs CV errors and files for plotting each K in r \
Use ```admixture_plot.R``` to plot figure in R 

**IQtree and RAXML**

IQtree and RAXML need to be installed in conda\
Run\
```iqtree.sh``` \
Run\
```raxml.sh``` 

**SVDQuartets**

generate file for SVDQuartets \
``` cat inputs/filtered_final_83_LD_pruned.recode.min4.nexus inputs/taxpartitions.txt > inputs/svd.nexus```\
PAUP needs to be installed and launched \
Run \
```exe inputs/svd.nexus``` \
Then run \
```svdq taxpartition=fish showScores=no seed=1234568 bootstrap nreps=1000 treeFile=svd.tre;```
To save the consenses tree in newick format run \
```savetree```

**SNAPP**

Beast2 v2.7.7 needs to be downloaded \
Run 2 times\
```scripts/SNAPP.sh``` 

**BFD***
Beast2 v2.6.7 needs to be downloaded \
Run\
```scripts/BFD.sh``` 

**DELINEATE**

DELINEATE needs to be downloaded to conda\
Guide tree is generated from the two SNAPP runs consensus tree using the TreeAnnotator application\
Run ```scripts/delineate.sh```

**TreeMix**

compress vcf file \
```bgzip inputs/filtered_final_94.recode.vcf > inputs/filtered_final_94.recode.vcf.gz``` \
index the compressed file \
```tabix inputs/filtered_final_94.recode.vcf.gz``` \
Run python file ```ipyrad_converter_treemix.py``` \
Followed by ```ipyrad_treemix.py```
Save output in folder names "tree" within the data folder for further analysis in R\
To plot use python file ```ipyrad_treemix_plot.py```

**Dsuite**

Download dsuite https://github.com/millanek/Dsuite \
Uses svd.newick.txt output from SVDquartets \
Run ```scripts/dsuite.sh```
Outputs dsuite.fb.txt \
To plot use dtools.py, which is part of the dsuite software\
run\
```python3 dtools.py dsuite.fb.txt inputs/svd.newick.txt```

**Nucleotide Diversity**

ANGSD needs to be downloaded in conda \
Run \
```diversity.sh```
Plotting is done in R using ```Jenynsia_genetic_diversity.R```













