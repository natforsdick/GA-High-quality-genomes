#!/bin/bash -e

#SBATCH -J psmc
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G 
#SBATCH -t 00:45:00 
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# converting consensus to psmc format, and testing parameter space

OUTDIR=/PATH/TO/rata-pop-gen/data/psmc/
cd $OUTDIR

module purge; module load psmc/0.6.5-gimkl-2018b

for fastq in EXT049-08_S8_map_rata_tahae_ref_diploid.fq.gz
do
    filename=$(basename "$fastq") 
    filename=${filename%.fq.gz}
    # Let's convert the diploid genome to PSMC suitable format
    if [ ! -e  ${filename}.psmcfa ]; then  
        echo "converting $filename"  
        fq2psmcfa -q20 ${fastq} > ${filename}.psmcfa
    else
        echo "$filename psmcfa file exists"
    fi

    # Running PSMC with standard parameters to get initial idea about how the parameters are working.
    echo "running psmc for $filename"
# -N max iterations, -t max 2N0 coalescent time, r initial theta/rho ratio - based on those used in Nadachowska-Brzyska et al 2015
    psmc -N50 -t11 -r5 -p "4+25*2+4+6" -o ${filename}-v4.psmc ${filename}.psmcfa
done
