#!/bin/bash

#SBATCH -J porechop
#SBATCH --time=01:20:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=58G
#SBATCH --output %x.%j.out 
#SBATCH --error %x.%j.err

# Adapter trimming/filtering of Nanopore data using porechop

# ENVIRONMENT #
module purge
module load Porechop/0.2.4-gimkl-2020a-Python-3.8.2

# Adapter removal
# If you intend to use Nanopolish downstream, you must use `--discard_middle` - this removes reads with adapters within them.
echo "Starting Porechop"
date

porechop -i /PATH/TO/INPUT/output/rata-MinION/Rata_2/sup-fastq/combined-sup-fastqs/rata-PB5-pass.fastq.gz \
    -o /PATH/TO/OUTPUT/output/rata-MinION/Rata_2/sup-fastq/combined-sup-fastqs/02-trimfilt/rata-adaprem.fastq \
    --discard_middle -t $SLURM_CPUS_PER_TASK
