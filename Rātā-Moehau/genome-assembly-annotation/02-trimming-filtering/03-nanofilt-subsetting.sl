#!/bin/bash

#SBATCH -J nanofilt-subset
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --output %x.%j.out 
#SBATCH --error %x.%j.err

# Remove short reads - can try several minimum lengths here

###############
# ENVIRONMENT #
module purge
module load nanofilt/2.6.0-gimkl-2020a-Python-3.8.2

cd /PATH/TO/INPUT/output/rata-MinION/Rata_2/sup-fastq/combined-sup-fastqs/02-trimfilt-b/

gunzip rata-trimmed.fastq.gz

# set minimum read lengths
subset="5000 10000"
for sub in $subset;
do
  echo "Starting NanoFilt for" $sub
  NanoFilt -l ${sub} rata-trimmed.fastq | gzip > rata-trimmed-${sub}.fastq.gz
  zcat rata-trimmed-${sub}.fastq.gz | echo $((`wc -l`/4))
done
