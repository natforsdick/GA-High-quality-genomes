#!/bin/bash -e

#SBATCH --job-name medaka
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

##########
# PARAMS #
ONTDIR=/PATH/TO/ONT/output/rata-MinION/
REFDIR=/PATH/TO/RACON-OUT/output/polish/racon/
REF=01-asm-shasta-purged-racon2.fasta
ONT=rata-all-trimmed.fastq
OUTDIR=/PATH/TO/OUTPUTS/output/polish/medaka/
##########

mkdir -p $OUTDIR && cd $OUTDIR
ml purge && ml medaka/1.11.1-Miniconda3-22.11.1-1

# we used r9.4.1 flowcells, on the MinION, with the sup model in Guppy v 6.2.1, so the closest model available is: r941_min_sup_g507
medaka_consensus -i ${ONTDIR}${ONT} -d ${REFDIR}${REF} -o ${OUTDIR} -t $SLURM_CPUS_PER_TASK -m r941_min_sup_g507
