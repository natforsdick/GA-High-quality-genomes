#!/bin/bash -e

# classifying unknowns (-u) in iterative fashion: run with 3 threads/cores (-t) and using only the known elements (-k) from the 
# same reference genome; append newly identified elements to the existing known element library (-a) and 
# write results to an output directory (-o). No Repbase classification is used here.
# We iterate over this several times, and stop when we see no extra gain in classifying unknown repeats (here we demonstrate up to 9 iterations, but typically 4-6 is enough).

ml purge && ml RepeatMasker/4.1.0-gimkl-2020a SeqKit/2.4.0 bioawk/1.0 RMBlast/2.10.0-GCC-9.2.0

cd /PATH/TO/OUTPUTS/output/07-annotation/repeats/

for i in {1..9}
do
j=$(expr $i + 1)
echo iteration $i

/PATH/TO/repclassifier -t 3 -u round-${i}_RepbaseEukaryota-Self/round-${i}_RepbaseEukaryota-Self.unknown \
-k round-${i}_RepbaseEukaryota-Self/round-${i}_RepbaseEukaryota-Self.known \
-a round-${i}_RepbaseEukaryota-Self/round-${i}_RepbaseEukaryota-Self.known -o round-${j}_RepbaseEukaryota-Self
done
