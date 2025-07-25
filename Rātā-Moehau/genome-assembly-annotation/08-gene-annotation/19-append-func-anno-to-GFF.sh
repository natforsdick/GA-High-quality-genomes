#!/bin/bash -e

# Replace semicolons delimiting different terms/descriptions with pipes
# Fix attribute keys to lowercase (uppercase attributes are reserved for official use in GFFs)
# Conditionally print GO field only when GO terms exist
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 {
  gsub(/;/, ",", $2) # Replace ; in PFAM with ,
  gsub(/;/, "|", $4)  # Replace semicolons with pipes in Description
  gsub(/;/, "|", $5)  # Replace semicolons with pipes in GO terms
  printf "%s\tpfam=%s;uniprot=%s;description=%s", $1, $2, $3, $4
  if ($5 != "" && $5 != "-") {
    printf ";go=%s", $5
  }
  print ""
}' combined_annotations_final.tsv > annotations_lookup.tsv

# Append functional annotations to GFF based on matching transcript_id= or ID=
awk -F'\t' 'BEGIN {
  OFS="\t"
  while ((getline < "annotations_lookup.tsv") > 0) {
    annots[$1] = $2
  }
}
{
  if ($3 == "transcript" && match($9, /ID=([^;]+)/, m)) {
    id = m[1]
    if (id in annots) {
      $9 = $9 ";" annots[id]
    }
  }
  print
}' ../tsebra-long-iso-startstop.gff > metBart-annotated-transcripts.gff

# create sequence-region lines to append to the GFF header based on the assembly scaffold info
cut -f1,2 /PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/05_full_out/metBart-contam-excl.simple_mask.soft.complex_mask.soft.fasta.fai > scaffold_lengths.tsv
awk '{print "##sequence-region", $1, 1, $2}' scaffold_lengths.tsv > sequence_regions.txt

# prepend seq-region lines
{
  # Extract the first line (GFF version header)
  head -n 1 metBart-annotated-transcripts.gff

  # Then add sequence-region lines
  cat sequence_regions.txt

  # Then add the rest of the original GFF (excluding first line)
  tail -n +2 metBart-annotated-transcripts.gff
} > metBart-annotated-transcripts-seqregions.gff

# We also want to combine this with the repeat annotations to produce one GFF containing all annotation types

# First we need to modify the repeat masked GFF, as some lines don't have coordinates in the target attribute values
sed 's/Target=/repeat_family=/g' /PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/05_full_out/metBart-contam-excl.full_mask.gff3 > repeat_annotations.cleaned.gff

# Retain headers from gene annotation GFF, then add the gene annotation lines, then grab just the repeat annotation lines from the full repeat mask GFF
grep '^#' metBart-annotated-transcripts-seqregions.gff > metBart-full-annotation.gff
grep -v '^#' metBart-annotated-transcripts-seqregions.gff >> metBart-full-annotation.gff
grep -v '^#' repeat_annotations.cleaned.gff >> metBart-full-annotation.gff
# then we want to sort based on position of the feature in the genome, so repeats and genes are properly interspersed
(grep '^#' metBart-full-annotation.gff; grep -v '^#' metBart-full-annotation.gff | sort -k1,1 -k4,4n) > metBart-full-annotation-sorted.gff

# Finally, let's check that the output GFF containing all the annotations is in a valid GFF format (we used this throughtout the steps above to ensure smooth sailing)
ml purge; ml genometools/1.6.1-GCCcore-7.4.0
gt gff3validator metBart-full-annotation-sorted.gff
# input is valid GFF3
