#!/bin/bash -e

#SBATCH -J psmc-ref-prep
#SBATCH --cpus-per-task=6
#SBATCH --mem=6G
#SBATCH -t 15:00:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# getting the files in the right format and location for psmc
# have copied the ref genome and the BAM alignment file into a new directory: /PATH/TO/rata-pop-gen/data/psmc
# here we use the unmasked reference genome that had potential contaminants removed through FCS-GX

REFDIR=/PATH/TO/OUTPUTS/output/scaffolding/yahs/
REF=genome.nextpolish2-mapped.PT_JBAT-out-postsynteny2.contam-excl.fa
INDIR=/PATH/TO/rata-pop-gen/data/psmc/

cd $INDIR

ml purge; ml BCFtools/1.19-GCC-11.3.0

echo "creating consensus for rata psmc"
# upper and lower coverage thresholds are set with the -d and -D flags
bcftools mpileup --threads 12 -C50 -f ${REF} EXT049-08_S8.aligned.sorted.bam | bcftools call --threads 12 -c - | vcfutils.pl vcf2fq -d 13 -D 90 | gzip > EXT049-08_S8_map_rata_tahae_ref_diploid.fq.gz  
