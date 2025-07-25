#!/bin/bash -e

# first we do data wrangling for the TSEBRA merged annotation output

REFDIR=/PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/05_full_out/
REF=metBart-contam-excl.simple_mask.soft.complex_mask.soft.fasta

cd /PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/merged-annotation/

# convert merged output GTF to GFF
ml purge; ml AGAT/1.0.0-gimkl-2022a-Perl-5.34.1-R-4.2.1

agat_convert_sp_gxf2gxf.pl -g tsebra.gtf -o tsebra.gff

# get some stats about the merged annotation 
# -g = assembly length in bp
agat_sq_stat_basic.pl -i tsebra.gff -g 279317961

# convert GFF to protein FASTA so we can compare compleasm BUSCO presence
ml purge; ml gffread/0.12.7-GCC-11.3.0
gffread -g $REFDIR$REF -y tsebra-protein.faa tsebra.gff

# split the protein FASTA so we can run BLAST all vs all as an array for faster processing
ml purge; ml genometools/1.6.1-GCCcore-7.4.0
gt splitfasta -numfiles 10 tsebra-protein.faa

# make sample list for BLAST input
ls tsebra-protein.faa.* >> fasta.list
