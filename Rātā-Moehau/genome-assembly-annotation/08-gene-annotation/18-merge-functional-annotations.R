# Here we join the functional annotation data generated from DIAMOND and InterProScan into a single file which will be used to annotate the GFF

# libraries
library(dplyr)
library(tidyr)
library(readr)

# Read InterProScan annotations
ipscan <- read_tsv("interpro_grouped_dedup.tsv", col_names = c("Transcript_ID", "PFAM_ID", "Description_IPScan", "GO_IPScan"))

# Replace pipes with semicolons in InterProScan GO terms
ipscan$GO_IPScan <- gsub("\\|", ";", ipscan$GO_IPScan)

# Read DIAMOND annotations
diamond <- read_tsv("merged_annotations_cleaned.sorted.tsv", col_names = c("Transcript_ID", "UniProt_ID", "Best_Hit", "Description_DIAMOND", "GO_DIAMOND"))

# Clean Diamond GO terms to remove space separators
diamond <- diamond %>%
  mutate(GO_DIAMOND = gsub(";\\s+", ";", GO_DIAMOND)) %>%  # Remove space after ;
  mutate(GO_DIAMOND = trimws(GO_DIAMOND))                  # Trim leading/trailing spaces

# Full join on Transcript_ID
merged <- full_join(ipscan, diamond, by = "Transcript_ID")

# Merge description and GO columns
# Function to combine and deduplicate semicolon-separated terms
combine_unique <- function(x, y) {
  all_vals <- unique(na.omit(c(unlist(strsplit(x, ";")), unlist(strsplit(y, ";")))))
  paste(all_vals[all_vals != ""], collapse = ";")
}

# Apply to descriptions and GO terms
merged <- merged %>%
  rowwise() %>%
  mutate(
    Description = combine_unique(Description_IPScan, Description_DIAMOND),
    GO_terms = combine_unique(GO_IPScan, GO_DIAMOND)
  ) %>%
  ungroup()

head(merged)

# Select and reorder final columns

final_output <- merged %>%
  select(Transcript_ID, PFAM_ID, UniProt_ID, Description, GO_terms)

# Write to file

write_tsv(final_output, "combined_annotations_final.tsv")
