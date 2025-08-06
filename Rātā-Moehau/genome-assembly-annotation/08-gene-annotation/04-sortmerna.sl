#!/bin/bash -e

#SBATCH -J sortmerna
#SBATCH --cpus-per-task=10
#SBATCH --mem=22G
#SBATCH -t 02:30:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# screen RNAseq data for data originating from rRNAs

##########
# PARAMS
INDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/02-kraken2/
OUTDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/03-sortmerna-out/
REFDIR=/PATH/TO/OUTPUTS/output/07-annotation/nat-pipeline/repeats/masking/05_full_out/
REF=metBart-contam-excl.simple_mask.soft.complex_mask.soft.fasta
FASTQ="Rata01AG1183009-clean-seqs Rata02AG1183010-clean-seqs"
mkdir -p $OUTDIR
cd $OUTDIR

ml purge; ml SortMeRNA/4.3.6

for FILE in ${FASTQ}
do
echo assessing for $FILE
sortmerna --ref ${REFDIR}${REF} --fastx --paired-out --reads ${INDIR}${FILE}_1.fq --reads ${INDIR}${FILE}_2.fq --workdir $OUTDIR --threads 18
