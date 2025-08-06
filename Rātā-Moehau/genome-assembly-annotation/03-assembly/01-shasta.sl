#!/bin/bash -e

#SBATCH -J shasta
#SBATCH --time=03:00:00
#SBATCH -c 36
#SBATCH --mem=150G
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err

# Assembling processed ONT data with shasta

###########
# ENVIRONMENT
SHASTA=/PATH/TO/SHASTA/shasta/shasta-Linux-0.10.0
INDIR=/PATH/TO/INPUT/output/rata-MinION/
OUTDIR=/PATH/TO/ASM/output/03-assembly/asm-shasta
CONFIG=/PATH/TO/CONFIG/scripts/rata-moehau-genome-assembly/03-assembly/shasta.config

$SHASTA --input ${INDIR}rata-all-trimmed.fastq --assemblyDirectory ${OUTDIR} --config $CONFIG --command assemble --threads 32
