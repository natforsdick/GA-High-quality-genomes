#!/bin/bash -e

#SBATCH -J blast-ava
#SBATCH --cpus-per-task=12
#SBATCH --mem=1G
#SBATCH -t 00:10:00
#SBATCH --out %x.%j.%a.out
#SBATCH --err %x.%j.%a.err
#SBATCH --array 0-9

# running blast all vs all of the output protein fasta sequences
# this allows checking for any issues due to unmasked TEs etc
INDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/merged-annotation/

ml purge; ml BLAST/2.16.0-GCC-12.3.0

cd $INDIR
# -in input fasta of annotated sequences
# make the blast DB so you can chunk the inputs and run in parallel:
makeblastdb -in tsebra-protein.faa -dbtype prot -out tsebra-protein.faa_blast_db

# Do all vs all blast of protein sequences:
SAMPLE_LIST=($(<fasta.list))
SAMPLE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

blastp -db tsebra-protein.faa_blast_db -query ${SAMPLE} -num_threads 12 -outfmt 6 -out rata-blast-ava-${SLURM_ARRAY_TASK_ID}.tsv
