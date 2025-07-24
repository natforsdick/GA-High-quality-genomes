#!/bin/bash -e

#SBATCH --job-name=make-juice-in 
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:10:00 
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err 

# Generating .hic file from YaHS output for visualisation and manual curation in Juicebox

##########
# PARAMS 
YAHSJUICE=/PATH/TO/yahs/juicer
OUTDIR=/PATH/TO/OUTPUTS/output/scaffolding/yahs/
JUICER=/PATH/TO/juicer/scripts/juicer_tools.1.9.9_jcuda.0.8.jar
REF_DIR=/PATH/TO/OUTPUTS/output/scaffolding/
REF=genome.nextpolish2.fa
REFPRE=genome.nextpolish2
SCAF=genome.nextpolish2-mapped.PT

##########
cd $OUTDIR

echo generating contact map
$YAHSJUICE pre -a -o ${SCAF}_JBAT ${SCAF}.bin ${OUTDIR}${SCAF}_scaffolds_final.agp \
${REF_DIR}${REF}.fai > ${SCAF}_JBAT.log 2>&1
echo done step 1

echo running juicer_tools pre
java -jar -Xmx16G $JUICER pre ${SCAF}_JBAT.txt ${SCAF}_JBAT.hic.part \
	<(cat ${SCAF}_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')

mv ${SCAF}_JBAT.hic.part ${SCAF}_JBAT.hic
echo done step 2
