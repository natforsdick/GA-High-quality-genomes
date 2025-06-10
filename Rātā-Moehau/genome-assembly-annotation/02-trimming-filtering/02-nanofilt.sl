#!/bin/bash

#SBATCH -J nano-filt-trim
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output %x.%j.out 
#SBATCH --error %x.%j.err

# Quality trimming and filtering of Nanopore data

###############
# ENVIRONMENT #
module purge
module load nanofilt/2.6.0-gimkl-2020a-Python-3.8.2

INDIR=/PATH/TO/INPUT/output/rata-MinION/Rata_2/sup-fastq/combined-sup-fastqs/02-trimfilt/
OUTDIR=/PATH/TO/OUTPUT/output/rata-MinION/Rata_2/sup-fastq/combined-sup-fastqs/02-trimfilt-b/
mkdir $OUTDIR
cd $INDIR

# Trim and filter to remove poor quality ends and short reads
# l = minimum length, q = minimum average quality
zcat rata-adaprem.fastq.gz | NanoFilt -l 500 -q 10 --headcrop 20 --tailcrop 20 | gzip > ${OUTDIR}/rata-trimmed.fastq.gz

# count retained reads and bases
zcat ${OUTDIR}rata-trimmed.fastq.gz | echo $((`wc -l`/4))
zcat ${OUTDIR}rata-trimmed.fastq.gz | awk 'NR%4==2 {sum += length($0)} END {print sum}' 
