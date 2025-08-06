#!/bin/bash -e

#SBATCH -J merge-tsebra
#SBATCH --cpus-per-task=2
#SBATCH --mem=3G
#SBATCH -t 00:10:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# Running TSEBRA to merge the outputs of Braker RNAseq and ODBViridiplantae+Myrtaceae annotations

INDIR1=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/braker/rna/
INDIR2=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/braker/orthodb/
OUTDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/merged-annotation/

mkdir -p $OUTDIR
cd $OUTDIR

ml purge; ml TSEBRA/1.1.2.5-gimkl-2022a-Python-3.11.3

tsebra.py -g ${INDIR1}Augustus/augustus.hints.gtf,${INDIR1}GeneMark-ET/genemark.gtf,${INDIR2}Augustus/augustus.hints.gtf,${INDIR2}GeneMark-EP/genemark.gtf \
    -e ${INDIR1}hintsfile.gff,${INDIR2}hintsfile.gff -o tsebra.gtf 2> tsebra.log
