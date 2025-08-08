#!/bin/bash -e

# If samples were sequenced across multiple lanes, you will then want to merge the resulting BAMs
# do this before collecting stats
# you will need a file containing a list of the duplicate BAMs to start, one BAM per line

###########
# PARAMS

BAMDIR=/PATH/TO/rata-pop-gen/data/output/04-mapped/bam/

# MODULES
module purge; module load SAMtools/1.10-GCC-9.2.0

cd $BAMDIR

# first, merge sequencing duplicates from the two sequencing batches
echo merging sequencing duplicates
mkdir duplicates
for duplicate in ${cat duplicate-list.txt}
do
  cp $duplicate duplicates/
done

# Then to merge them, you will need to replace 'samp' below with the specific files for merging, where the *-merged.aligned.sorted.bam is the output file
# repeat as many times as required
cd duplicates/
samtools merge ../samp01-merged.aligned.sorted.bam samp01-S01.aligned.sorted.bam samp01-S02.aligned.sorted.bam
