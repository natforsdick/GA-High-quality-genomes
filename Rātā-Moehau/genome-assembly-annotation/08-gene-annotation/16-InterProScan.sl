#!/bin/bash -e

#SBATCH -J ipscan
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH --time=01:30:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# running InterProScan for functional annotation of merged processed annotation outputs

ml purge; ml InterProScan/5.66-98.0-gimkl-2022aPerl-5.34.1-Python-3.11.3

FASTA=metBart-contam-excl-prot_clean.fasta
OUTDIR=/PATH/TO/OUTPUTS/output/07-annotation/gene-annotation/merged-annotation/functional/
TMPDIR=/PATH/TO/tmp-IPScan/

mkdir -p $OUTDIR
mkdir -p $TMPDIR
cd $OUTDIR

# interproscan.sh -i [proteins.fasta] -t [p/n = prot or nucl input] -appl Pfam -cpu [cpu] -d [outdir] \
# -goterms -T [tmpdir] --pathways - gets pathway annotations
interproscan.sh -i ${INDIR}${FASTA} -t p -appl Pfam -cpu 10 -d $OUTDIR -goterms -T $TMPDIR --pathways
