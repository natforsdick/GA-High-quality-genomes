# combine the UniProt best hits information from DIAMOND with the output from the UniProt accession ID mapping

# first we need to sort by accession ID
awk '{
  split($2, a, "|");
  print $1 "\t" a[2];
}' metBart-contam-excl-prot-UniPSP_best-hits.tsv > best_hits_clean.tsv

sort -k2,2 best_hits_clean.tsv > best_hits.sorted.tsv

# Extract header from the ID mapping list (first line)
head -n 1 idmapping.tsv > idmapping_header.tsv

# Sort all lines except the header
tail -n +2 idmapping.tsv | sort -k1,1 > idmapping_sorted.tsv

# Recombine header and sorted data
cat idmapping_header.tsv idmapping_sorted.tsv > annotations_sorted.tsv
tail -n +2 annotations_sorted.tsv > annotations_sorted_noheader.tsv
join -1 2 -2 1 -t $'\t' \
  -o 1.1 1.2 2.2 2.5 2.10 \
  best_hits.sorted.tsv annotations_sorted_noheader.tsv > merged_annotations.tsv

# we don't actually need the info stored in brackets from the ID mapping (gene aliases), so let's remove this
# Remove all square or round brackets and their contents from description
awk -F'\t' -v OFS='\t' '{
    gsub(/\[[^]]*\]/, "", $4)  # Remove [ ... ]
    gsub(/\([^)]*\)/, "", $4)  # Remove ( ... )
    print
}' merged_annotations.tsv > merged_annotations_cleaned.tsv

# then we need to process the InterProScan outputs
