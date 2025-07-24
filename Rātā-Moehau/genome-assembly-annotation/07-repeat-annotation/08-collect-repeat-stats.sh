#!/bin/bash -e

# summarise repeat metrics

##########
# PARAMS
INDIR=PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/
OUTDIR=/PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/05_full_out/
REF=metBart-contam-excl

##########
mkdir -p $OUTDIR
cd $INDIR

# combine full RepeatMasker result files - .cat.gz
cat 01_simple_out/${REF}.simple_mask.cat.gz \
02_eukaryota_out/${REF}.eukaryota.masked.fasta.cat.gz \
03_known_out/${REF}.known.masked.cat.gz \
04_unknown_out/${REF}.unknown.masked.fasta.cat.gz \
> 05_full_out/${REF}.full_mask.cat.gz

# combine RepeatMasker tabular files for all repeats - .out
cat 01_simple_out/${REF}.simple_mask.out \
<(cat 02_eukaryota_out/${REF}.eukaryota.masked.fasta.out | tail -n +4) \
<(cat 03_known_out/${REF}.known.masked.out | tail -n +4) \
<(cat 04_unknown_out/${REF}.unknown.masked.fasta.out | tail -n +4) \
> 05_full_out/${REF}.full_mask.out

# copy RepeatMasker tabular files for simple repeats - .out
cat 01_simple_out/${REF}.simple_mask.out > 05_full_out/${REF}.simple_mask.out

# combine RepeatMasker tabular files for complex, interspersed repeats - .out
cat 02_eukaryota_out/${REF}.eukaryota.masked.fasta.out \
<(cat 03_known_out/${REF}.known.masked.out | tail -n +4) \
<(cat 04_unknown_out/${REF}.unknown.masked.fasta.out | tail -n +4) \
> 05_full_out/${REF}.complex_mask.out

# combine RepeatMasker repeat alignments for all repeats - .align
cat 01_simple_out/${REF}.simple_mask.align \
02_eukaryota_out/${REF}.eukaryota.masked.fasta.align \
03_known_out/${REF}.known.masked.align \
04_unknown_out/${REF}.unknown.masked.fasta.align \
> 05_full_out/${REF}.full_mask.align

##########
# resummarize repeat compositions from combined analysis of all RepeatMasker rounds
ml purge && ml RepeatModeler/2.0.3-Miniconda3
ProcessRepeats -a -species eukaryota 05_full_out/${REF}.full_mask.cat.gz 2>&1 | tee logs/05_fullmask.log

##########
# PARAMS
REFDIR=/PATH/TO/OUTPUTS/output/07-annotation/
REF=metBart-contam-excl

# calculate the length of the genome sequence in the FASTA
ml purge && ml seqtk/1.4-GCC-11.3.0 Perl/5.34.1-GCC-11.3.0
seqtk comp ${REFDIR}${REF}.fasta > refstats.txt
allLen=`awk '{sum+=$2;} END{print sum;}' refstats.txt`
# calculate the length of the N sequence in the FASTA
nLen=`seqtk comp ${REFDIR}${REF}.fasta | datamash sum 9` 

# tabulate repeats per subfamily with total bp and proportion of genome masked
cat 05_full_out/${REF}.full_mask.out | tail -n +4 | awk -v OFS="\t" '{ print $6, $7, $11 }' |\
awk -F '[\t/]' -v OFS="\t" '{ if (NF == 3) print $3, "NA", $2 - $1 +1; else print $3, $4, $2 - $1 +1 }' |\
datamash -sg 1,2 sum 3 | grep -v "\?" |\
awk -v OFS="\t" -v genomeLen="${allLen}" '{ print $0, $3 / genomeLen }' > 05_full_out/${REF}.full_mask.tabulate
