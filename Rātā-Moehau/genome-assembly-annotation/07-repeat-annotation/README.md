# Repeat annotation

Once all QC has been completed for the assembly, including excluding of any contaminating sequences via FCS-GX, we are ready to perform annotation.

First, we annotate repeats via an iterative process. This process uses RepeatModeler and RepeatMasker, and follows the process described by Daren Card [here](https://darencard.net/blog/2022-07-09-genome-repeat-annotation/).
The end result will be three GFFs containing information on all repeats, only simple repeats, and only complex repeats, and a version of the assembly in FASTA format where all repeats are soft-masked, which will be used as input for gene annotation. 
