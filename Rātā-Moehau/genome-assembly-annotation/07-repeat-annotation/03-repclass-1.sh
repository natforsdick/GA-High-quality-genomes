#!/bin/bash -e

# Classifying repeatmodeler output unknowns (-u): run with 3 threads/cores (-t) and using the eukaryote elements (-d) from Repbase
# and known elements (-k) from the same reference genome.
# Append newly identified elements to the existing known element library (-a) and write results to an output directory (-o).

ml purge && ml RepeatMasker/4.1.0-gimkl-2020a SeqKit/2.4.0 bioawk/1.0 RMBlast/2.10.0-GCC-9.2.0

TMPDIR=/nesi/nobackup/landcare03691/tmp-repclass/
mkdir -p $TMPDIR
export TMPDIR=$TMPDIR

cd /PATH/TO/OUTPUTS/output/07-annotation/repeats/
/PATH/TO/repclassifier -t 3 -d eukaryota -u metBart-families.metBart1.fa.unknown \
  -k metBart-families.metBart1.fa.known -a metBart-families.metBart1.fa.known -o round-1_RepbaseEukaryota-Self
