#!/bin/bash -e

#SBATCH -J repmask3
#SBATCH --cpus-per-task=20
#SBATCH --mem=18G
#SBATCH -t 1:30:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# rounds 3-4 of RepeatMasker to annotate known and unknown elements from the species-specific de novo repeat library

##########
# PARAMS
OUTDIR=/PATH/TO/OUTPUTS/output/07-annotation/repeats/masking/
REFDIR=/PATH/TO/OUTPUTS/output/07-annotation/
REF=metBart-contam-excl

##########
ml purge && ml RepeatMasker/4.1.0-gimkl-2020a

cd $OUTDIR

##########
# round 3: annotate/mask known elements sourced from species-specific de novo repeat library using output froom 2nd round of RepeatMasker
echo starting round 3
RepeatMasker -pa 24 -a -e ncbi -dir 03_known_out -nolow \
-lib ../round-10_RepbaseEukaryota-Self/round-10_RepbaseEukaryota-Self.known \
02_eukaryota_out/${REF}.eukaryota.masked.fasta 2>&1 | tee logs/03_knownmask.log

# rename outputs
rename eukaryota known 03_known_out/${REF}*
rename .masked .fasta known 03_known_out/${REF}*

##########
# round 4: annotate/mask unknown elements sourced from species-specific de novo repeat library using output froom 3nd round of RepeatMasker
echo starting round 4
RepeatMasker -pa 24 -a -e ncbi -dir 04_unknown_out -nolow \
-lib ../round-10_RepbaseEukaryota-Self/round-10_RepbaseEukaryota-Self.unknown \
03_known_out/${REF}.known.masked.fasta 2>&1 | tee logs/04_unknownmask.log

# rename outputs
rename .known .unknown 04_unknown_out/${REF}*
rename .masked .fasta known 04_unknown_out/${REF}*
