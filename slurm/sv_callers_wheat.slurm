#!/bin/bash
#SBATCH --job-name=sv_callers_wheat
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/sv_calling_moryzae/wheat_blast_11_2_22

SAMPLE_LIST=wheat_blast_illumina_seq_runs.txt
REFERENCE=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_to_annotate/GCA_004785725.2_ASM478572v2_genomic.fna


if [ ! -f "${REFERENCE}.amb" ]; then
    bwa index ${REFERENCE}
fi

if [ -f slurm_ids ]; then
    rm slurm_ids
fi

while read SAMPLE; do
    if [ ! -d "$SAMPLE" ]; then
        mkdir $SAMPLE
    fi
    if [ ! -d "${SAMPLE}/data/" ]; then
        mkdir ${SAMPLE}/data/
    fi
    if [ ! -s "${SAMPLE}/data/${SAMPLE}.bam" ]; then
        echo "starting bwa mem for ${SAMPLE}"
        cd ${SAMPLE}/data/
            sbatch --parsable --job-name=${SAMPLE}.sra_map --export=SAMPLE=$SAMPLE,REFERENCE=$REFERENCE /global/home/users/pierrj/git/slurm/download_sra_and_map.slurm >> ../../slurm_ids
        cd ../..
    fi
done < $SAMPLE_LIST

if [ -f slurm_ids ]; then
    # wait for all jobs to finish
    sbatch -W --dependency=afterok:$(cat slurm_ids | head -c -1 | tr '\n' ':') /global/home/users/pierrj/git/slurm/dummy_job.slurm

    wait
fi

# run all callers with per sample parallelization

module purge

module load bcftools

for caller in wham lumpy manta delly; do
    if [ -f ${caller}_jobqueue ]; then
        rm ${caller}_jobqueue
    fi

    while read SAMPLE; do
        if [ ! -d "${SAMPLE}/${caller}/" ]; then
            mkdir ${SAMPLE}/${caller}/
        fi
        if [ ! -s "${SAMPLE}/${caller}/${SAMPLE}.${caller}.vcf" ]; then
            echo /global/home/users/pierrj/git/bash/run_${caller}.sh -f ${SAMPLE}/data/${SAMPLE}.bam -r ${REFERENCE} -t 1 -o ${SAMPLE}/${caller}/${SAMPLE}.${caller}.vcf -d ${SAMPLE}/${caller}/ >> ${caller}_jobqueue
        fi
        
    done < $SAMPLE_LIST

    if [[ ! "${caller}" == "wham" ]]; then
        source activate /global/scratch/users/pierrj/conda_envs/${caller}
        export PERL5LIB=/global/scratch/users/pierrj/conda_envs/${caller}/lib
    fi

    if [ -f ${caller}_jobqueue ]; then
        echo "starting ${caller} jobs"
        parallel -j ${SLURM_NTASKS} < ${caller}_jobqueue
    fi

done

## merge all

## in order
# distance between breakpoints
# number of supporting callers
# take type into account
# take strand into account
# disabled?
# minimum size
# output

while read SAMPLE; do
    if [ ! -d "${SAMPLE}/all/" ]; then
        mkdir ${SAMPLE}/all/
    fi

    if [ -f  ${SAMPLE}/all/${SAMPLE}.survivor.vcflist ]; then
        rm  ${SAMPLE}/all/${SAMPLE}.survivor.vcflist
    fi

    echo ${SAMPLE}/delly/${SAMPLE}.delly.vcf >> ${SAMPLE}/all/${SAMPLE}.survivor.vcflist
    echo ${SAMPLE}/lumpy/${SAMPLE}.lumpy.vcf >> ${SAMPLE}/all/${SAMPLE}.survivor.vcflist
    echo ${SAMPLE}/manta/${SAMPLE}.manta.vcf >> ${SAMPLE}/all/${SAMPLE}.survivor.vcflist
    echo ${SAMPLE}/wham/${SAMPLE}.wham.vcf >> ${SAMPLE}/all/${SAMPLE}.survivor.vcflist

    /global/scratch/users/pierrj/SURVIVOR/Debug/SURVIVOR merge ${SAMPLE}/all/${SAMPLE}.survivor.vcflist 1000 3 1 1 0 50 ${SAMPLE}/all/${SAMPLE}.all.vcf

done < $SAMPLE_LIST