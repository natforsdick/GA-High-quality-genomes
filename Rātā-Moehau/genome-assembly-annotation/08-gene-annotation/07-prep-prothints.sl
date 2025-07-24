#!/bin/bash -e

#SBATCH -J prep-prot
#SBATCH -c 18
#SBATCH --mem=18G
#SBATCH --time=08:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# Preparing the OrthoDB Viridiplantae + genbank myrtaceae protein file for input to BRAKER
# Here we pulled down the protein.faa's for those Myrtaceae with annotation information available, 
# and appended that to the Viridiplantae OrthoDB

# First, download protein.faa files - downloaded 2025/09/05.
# Then we concatenate these 
cd /PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/reference-annotations/
cat GC*.faa > genbank-20240905-myrtaceae.faa
cat Viridiplantae.fa genbank-20240905-myrtaceae.faa > ODBViridiplantae-genbank-20240905-myrtaceae.faa

# Then we need to remove sequences that contain any characters that BRAKER can't handle
grep "[?*<&%#@]" -B1 ODBViridiplantae-genbank-20240905-myrtaceae.faa > badchar.list
grep "^>" badchar.list > badchar.seq.list
wc -l badchar.seq.list # check how many sequences this is

# manually remove those sequences and check that output is correctly formatted
ml purge; ml bioawk
bioawk -cfastx '{printf(">%s\t%s\n", $name, $seq)}' ODBViridiplantae-genbank-20240905-myrtaceae.faa |\
grep -v -f badchar.seqhead.list |\
tr "\t" "\n" > ODBViridiplantae-genbank-20240905-myrtaceae.clean.faa

# now we are ready to process this input for BRAKER
ml purge; ml BRAKER/3.0.8-gimkl-2022a-Perl-5.34.1

prothint.py /PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/05_full_out/metBart-contam-excl.simple_mask.soft.complex_mask.soft.fasta \
  /PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/reference-annotations/ODBViridiplantae-genbank-20240905-myrtaceae.clean.faa --threads 32
