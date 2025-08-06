#!/bin/bash -e

#SBATCH -J repmask1
#SBATCH --cpus-per-task=26
#SBATCH --mem=18G
#SBATCH -t 01:00:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# RepeatMasker round 1: annotate/mask simple repeats

##########
# PARAMS
OUTDIR=/PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/
REFDIR=/PATH/TO/OUTPUTS/output/07-annotation/
REF=metBart-contam-excl.fasta
base=$(basename "$REF")
base=${base%.*}

##########
ml purge && ml RepeatMasker/4.1.0-gimkl-2020a

mkdir -p ${OUTDIR}logs
cd $OUTDIR

RepeatMasker -pa 24 -a -e ncbi -dir ${OUTDIR}01_simple_out/ -noint -xsmall ${REFDIR}${REF} 2>&1 | tee ${OUTDIR}logs/01_simplemask.log

# renaming outputs
rename fasta simple_mask ${OUTDIR}01_simple_out/${base}*
rename .masked .masked.fasta 01_simple_out/${base}*
