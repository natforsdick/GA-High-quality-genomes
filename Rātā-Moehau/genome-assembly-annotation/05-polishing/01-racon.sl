#!/bin/bash -e

#SBATCH --job-name racon
#SBATCH --mem=105G
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition gpu,hgx
#SBATCH --gpus-per-node A100:1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

##########
# PARAMS #
ONTDIR=/PATH/TO/PROCESSED-ONT/output/rata-MinION/
REFDIR=/PATH/TO/PURGED-ASM/output/04-purge-dups/asm-shasta/
REF=01-asm-shasta-purged
ONT=rata-all-trimmed.fastq
OUTDIR=/PATH/TO/OUTPUTS/output/polish/racon/
##########

mkdir -p $OUTDIR && cd $OUTDIR

module purge && module load {CUDA/11.4.1,Racon/1.5.0-GCC-11.3.0,minimap2/2.24-GCC-11.3.0}
ml list

echo 1X Racon
#minimap2 -t 10 ${REFDIR}${REF}.fa ${ONTDIR}${ONT} > ${REF}-ONT.paf
racon -t 10 ${ONTDIR}${ONT} ${REF}-ONT.paf ${REFDIR}${REF}.fa > ${REF}-racon1.fasta

echo 2X Racon
minimap2 -t 10 ${REF}-racon1.fasta ${ONTDIR}${ONT} > ${REF}-racon1-ONT.paf
racon -t 10 ${ONTDIR}${ONT} ${REF}-racon1-ONT.paf ${REF}-racon1.fasta > ${REF}-racon2.fasta
echo 2x Racon completed
