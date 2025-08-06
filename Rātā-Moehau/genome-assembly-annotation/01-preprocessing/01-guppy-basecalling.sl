#!/bin/bash -e

#SBATCH --job-name       ont-rata-guppy
#SBATCH --gpus-per-node  A100:1
#SBATCH --mem            4G
#SBATCH --cpus-per-task  4
#SBATCH --time           11:00:00 
#SBATCH --output         %x.%j.out
#SBATCH --error         %x.%j.err

# Running Guppy basecalling for Nanopore fast5 reads from ONT R9.4.1 flowcells.

#########
# PARAMS
INDIR=/PATH/TO/INPUT/fast5/
OUTDIR=/PATH/TO/output/rata-MinION/Rata_2/
SUPCFG=/PATH/TO/ONT-scripts/guppy-cfg/dna_r9.4.1_450bps_sup.cfg # set path to correct ONT config file corresponding to flowcell
#########

module purge
module load ont-guppy-gpu/6.2.1

guppy_basecaller -i ${INDIR} -s ${OUTDIR}sup-fastq --config ${SUPCFG} \
  --device auto --recursive --records_per_fastq 4000 \
  --detect_mid_strand_adapter
