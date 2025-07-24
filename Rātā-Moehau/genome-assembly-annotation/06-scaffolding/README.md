# Assembly scaffolding

Steps 1-2 are to QC the preliminary Omni-C data produced from a MiSeq run to test the quality of the library ahead of full sequencing. 

Steps 3 onwards are for the full sequencing run, if the preliminary data passes QC. 

In step 3, short-read data are cleaned with fastp.

Step 4 prepares the genome index files.

Step 5 follows the Dovetail™ Omni-C™ pipeline for aligning the Omni-C data to the assembly and processing the alignments. This pipeline implements BWA, SAMtools, and Pairtools. 
Following scaffolding and manual curation, a second round of aligment may be performed, which would require preprocessing as in step 4, and modification of the assembly inputs for step 5.

Step 6 sorts the alignments, and is presented as a standalone step due to the relatively higher memory requirements. 

Scaffolding is performed with YaHS in step 7.

Step 8 contains post-processing steps necessary for visualisation and manual curation of the scaffolded assembly in Juicebox. 

Following any manual curation, step 9 implements conversion from the Juicebox output to FASTA for final QC, and any subsequent repeat of the previous steps. 
