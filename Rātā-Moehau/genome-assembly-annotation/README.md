# Genome assembly and annotation

This directory contains all scripts required for assembly and annotation of the rātā Moehau genome, from pre-processing raw short- and long-read data, initial assembly, scaffolding, annotation, QC, and post-processing steps. 

We recommend that various tools are used for assessing the quality of the assembly at various steps in the process, including gathering standard assembly metrics (assembly length, contig/scaffold length, contig/scaffold N50, distribution of contig/scaffold lengths, number/length of gaps), ortholog presence via BUSCO or Compleasm, alignment of short-/long-read data to assess coverage spikes or gaps, and assessing assembly completeness with Merqury. We also used TIDK to look for potential telomeres, and FCS-GX to check for contamination. 
