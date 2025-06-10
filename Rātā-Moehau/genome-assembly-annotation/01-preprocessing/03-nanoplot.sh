#!/bin/bash -e

# Assessing sequence data with NanoPlot - need to install: <https://github.com/wdecoster/NanoPlot>
# Takes 1 param: input directory
# -p = output prefix, -o = output directory, -c = colour scheme, --summary = sequencing summary txt file produced by guppy

INDIR=/PATH/TO/INPUT/output/rata-MinION/Rata_2/sup-fastq/

# To make conda work in your shell env:
source /PATH/TO/Miniconda3/4.8.3/etc/profile.d/conda.sh

conda activate NanoPlot
cd $INDIR
NanoPlot -o /PATH/TO/OUTPUT/output/rata-MinION/Rata_2/sup-fastq/sup-QC/ -p nanoplot-rata-batch2 -c forestgreen --N50 --summary sequencing_summary.txt
conda deactivate
