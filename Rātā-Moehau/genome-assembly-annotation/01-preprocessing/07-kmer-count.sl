#!/bin/bash -e

#SBATCH -J jellyfish
#SBATCH --time=02:00:00 # 35 min is sufficient for jellyfish steps
#SBATCH -c 36
#SBATCH --mem=30G
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err

# Counting kmers for genome size estimation.
# Can use GenomeScope web assessment to estimate genome size.
# Requires jellyfish kmer histogram report as input.
# https://genome.umd.edu/docs/JellyfishUserGuide.pdf

##########
# PARAMS
DATADIR=/PATH/TO/INPUT/output/illumina/02-trimmed/
READS=EXT041-06_S6
mer=21
CPU=20
##########

ml purge
ml Jellyfish/2.3.0-gimkl-2020a
cd ${DATADIR}

# To generate: -C = canonical (necessary), -m = kmers, -s = hash size, -t = threads

echo 'Beginning jellyfish bloom'
date
jellyfish bc -C -m ${mer} -s 3G -t ${CPU} -o rata-${mer}mer.bc <(zcat ${READS}_L001_val_1.fq.gz) <(zcat ${READS}_L001_val_2.fq.gz) <(zcat ${READS}_L002_val_1.fq.gz) <(zcat ${READS}_L002_val_2.fq.gz)

echo 'Beginning jellyfish count'
jellyfish count -C -m ${mer} -s 3G -t ${CPU} --bc rata-${mer}mer.bc -o rata-${mer}mer-counts.jf <(zcat ${READS}_L001_val_1.fq.gz) <(zcat ${READS}_L001_val_2.fq.gz) <(zcat ${READS}_L002_val_1.fq.gz) <(zcat ${READS}_L002_val_2.fq.gz)

echo 'Beginning jellyfish histo'
date
jellyfish histo -t ${CPU} rata-${mer}mer-counts.jf > rata-${mer}mer-histo.out

echo 'Completed kmer-count'
date
