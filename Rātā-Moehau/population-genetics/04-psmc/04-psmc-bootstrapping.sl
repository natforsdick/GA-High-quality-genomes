#!/bin/bash -e

#SBATCH -J psmc-boot
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH -t 12:00:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# run bootstrapping for psmc

OUTDIR=/PATH/TO/rata-pop-gen/data/psmc/
cd $OUTDIR

module purge; module load psmc/0.6.5-gimkl-2018b

fastq=EXT049-08_S8_map_rata_tahae_ref_diploid.fq.gz
filename=$(basename "$fastq")
filename=${filename%.fq.gz}

# split to run smaller chunks
echo splitting $fastq
splitfa ${filename}.psmcfa > ${filename}_split.psmcfa

# perform bootstrapping - check that the parameter space is correct
echo bootstrapping $fastq
seq 100 | xargs -i -P 8 psmc -N50 -t11 -r5 -b -p "4+25*2+4+6" -o ${filename}_round-{}.psmc ${filename}_split.psmcfa | sh
