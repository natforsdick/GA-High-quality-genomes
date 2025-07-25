#!/bin/bash -e

# after looking at the metrics and some QC, especially the compleasm results, it's apparent that there is a large proportion of 
# duplication within the merged annotation.
# This likely results from the way TSEBRA merges the annotations. But we can do some processing to clean these up and improve the
# annotation quality at this stage.

REFDIR=/PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/05_full_out/
REF=metBart-contam-excl.simple_mask.soft.complex_mask.soft.fasta
ASM=$REFDIR$REF
LEN=279317961 # genome asm length in bp

cd /PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/merged-annotation/
ml purge; ml AGAT/1.0.0-gimkl-2022a-Perl-5.34.1-R-4.2.1

# first let's remove overlapping/nested genes from the annotation
agat_sp_fix_overlaping_genes.pl --gff tsebra.gff -o tsebra-rem-overlaps.gff  &> rem-overlaps.log

# now let's retain only the longest isoform for each transcript
agat_sp_keep_longest_isoform.pl -gff tsebra-rem-overlaps.gff -o tsebra-long-iso.gff &> longest-iso.log

# AGAT requires the reference assembly to be line wrapped in a particular way for this next step
ml purge; ml FASTX-Toolkit/0.0.14-GCC-11.3.0
fasta_formatter -i $ASM -o fold-metBart-contam-excl.simple_mask.soft.complex_mask.soft.fasta -w 80

# remove transcripts that are incomplete (missing start/stop codons)
ml purge; ml AGAT/1.0.0-gimkl-2022a-Perl-5.34.1-R-4.2.1
agat_sp_filter_incomplete_gene_coding_models.pl --gff tsebra-long-iso.gff \
  --fasta fold-metBart-contam-excl.simple_mask.soft.complex_mask.soft.fasta \
  -o tsebra-long-iso-startstop.gff &> startstop.log

# again, we want to gather some metrics at this stage, and assess this cleaned annotation with compleasm to see if this has improved things
agat_sq_stat_basic.pl -i tsebra-long-iso-startstop.gff -g $LEN -o tsebra-long-iso-startstop.basic
# -d = plot distribution
agat_sp_statistics.pl --gff tsebra-long-iso-startstop.gff --gs $LEN -d -o tsebra-long-iso-startstop.stat

# convert output gff to fasta for BUSCO screening
agat_sp_extract_sequences.pl -g tsebra-long-iso-startstop.gff -f fold-metBart-contam-excl.simple_mask.soft.complex_mask.soft.fasta \
  -p -o metBart-contam-excl-prot.fasta

# the output protein fasta includes '*' used to denote stop codons, but InterProScan can't read that format, so let's strip those out
sed 's/\*//g' metBart-contam-excl-prot.fasta > metBart-contam-excl-prot_clean.fasta
