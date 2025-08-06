#!/bin/bash -e
#SBATCH --job-name=diamond
#SBATCH --account=landcare03691
#SBATCH --time=00:05:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=3
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# Using protein blast in DIAMOND for functional annotation of processed protein sequences 
# Before running this script, be sure to download and unzip the latest UniProt/SwissProt database from:
# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

OUTDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/merged-annotation/functional/
PROTEIN=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/merged-annotation/metBart-contam-excl-prot.fasta # input protein file from annotation
UNIPROT=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/reference-annotations/uniprot_sprot.fasta
TMPDIR="/PATH/TO/tmp-diamond-${SLURM_JOB_ID}"

mkdir -p $TMPDIR
cd $OUTDIR
module purge; module load DIAMOND/2.1.10-GCC-12.3.0

# create a diamond-formatted version of the database - only need to do this once
diamond makedb --in ${UNIPROT} -d uniprot_sprot --threads 2

# run a search in blastp mode
diamond blastp -d uniprot_sprot --query ${PROTEIN} --out metBart-contam-excl-prot-UniPSP-matches-inf.tsv --threads 4 --tmpdir $TMPDIR --sensitive --evalue 0.00001 \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle
