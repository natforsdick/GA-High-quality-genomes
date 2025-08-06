#!/bin/bash -e

# We then want to do a little processing of the DIAMOND blast outputs in preparation for merging with InterProScan outputs

OUTDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/merged-annotation/functional/
cd $OUTDIR

# retain just the best hit for each transcript
sort -k1,1 -k12,12nr ${outdir}metBart-contam-excl-prot-UniPSP-matches-inf.tsv | awk '!seen[$1]++' > ${outdir}metBart-contam-excl-prot-UniPSP_best-hits.tsv

# check number retained matches the number of queries aligned by DIAMOND
wc -l metBart-contam-excl-prot-UniPSP_best-hits.tsv 

# Let's then extract the list of UniProt accession IDs so we can collect the functional annotation information
cut -f2 metBart-contam-excl-prot-UniPSP_best-hits.tsv | cut -d'|' -f2 > accessions.txt
# then map accession IDs to functional annotations at https://www.uniprot.org/id-mapping
