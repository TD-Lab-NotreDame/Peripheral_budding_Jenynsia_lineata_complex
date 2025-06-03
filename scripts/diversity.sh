#!/bin/bash
#$ -q long
#$ -pe smp 1
#$ -N diversity

echo "start"

# Reference genome fasta from ddocent output of 83 individuals dataset
REF="reference.fasta"

# Array of population names (matching your bamlist filenames)
pops=(lineataN lineataW luxata onca lineataE lienataEW darwini)

for pop in "${pops[@]}"
do
    echo "Processing population: $pop"

angsd -b inputs/${pop}.txt \
      -ref $REF \
      -anc $REF \
      -out ${pop}_angsd \
      -uniqueOnly 1 \
      -remove_bads 1 \
      -only_proper_pairs 1 \
      -minMapQ 30 \
      -minQ 20 \
      -baq 1 \
      -setMinDepth 10 \
      -setMaxDepth 58 \
      -GL 1 \
      -doCounts 1 \
      -doSaf 1

    # Calculate SFS from the SAF file
    realSFS ${pop}_angsd.saf.idx -fold 1 > ${pop}_folded.sfs

    # Calculate diversity statistics (theta) using SFS
    realSFS saf2theta ${pop}_angsd.saf.idx -sfs ${pop}_folded.sfs -fold 1 -outname ${pop}_angsd

    # Get sliding window estimates (adjust window size as needed)
    thetaStat do_stat ${pop}_angsd.thetas.idx
    thetaStat print ${pop}_angsd.thetas.idx > ${pop}_thetas_per_site.txt

done


echo "end"
