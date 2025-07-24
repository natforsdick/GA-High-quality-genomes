#!/bin/bash -e

#SBATCH -J repmask2
#SBATCH --cpus-per-task=24
#SBATCH --mem=10G
#SBATCH -t 1:10:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# round 2: annotate/mask complex, interspersed repeats

##########
# PARAMS
OUTDIR=/PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/
REFDIR=/PATH/TO/OUTPUTS/output/07-annotation/
REF=metBart-contam-excl.simple_mask.masked.fasta 

##########
ml purge && ml RepeatMasker/4.1.0-gimkl-2020a

cd $OUTDIR

RepeatMasker -pa 24 -a -e ncbi -dir ${OUTDIR}02_eukaryota_out/ -nolow \
-species eukaryota ${OUTDIR}01_simple_out/${REF} 2>&1 | tee ${OUTDIR}logs/01_eukaryota.log

# rename outputs
rename .simple_mask.masked.fasta .eukaryota.masked.fasta ${OUTDIR}02_eukaryota_out/metBart-contam-excl*
