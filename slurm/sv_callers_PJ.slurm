#!/bin/bash
#SBATCH --job-name=sv_callers_PJ
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


cd /global/scratch/users/pierrj/sv_calling_moryzae/run_1_5_2022

SAMPLE_LIST=test_sample_list
REFERENCE=/global/scratch/users/pierrj/references/moryzae_70-15_ref_with_mito.fasta # must be full path

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
    if [ ! -f "${SAMPLE}/data/${SAMPLE}.bam" ]; then
        cd ${SAMPLE}/data/
            sbatch --parsable --job-name=${SAMPLE}.sra_map --export=SAMPLE=$SAMPLE,REFERENCE=$REFERENCE /global/home/users/pierrj/git/slurm/download_sra_and_map.slurm >> slurm_ids
        cd ..
    fi
done < SAMPLE_LIST

# wait for all jobs to finish
sbatch -W --dependency=afterok:$(cat slurm_ids | head -c -1 | tr '\n' ':') /global/home/users/pierrj/git/slurm/dummy_job.slurm
wait

# run all callers with per sample parallelization

module load bcftools

for caller in wham lumpy manta delly; do
    if [ -f ${caller}_jobqueue ]; then
        rm ${caller}_jobqueue
    fi

    while read SAMPLE; do
        if [ ! -d "${SAMPLE}/${caller}/" ]; then
            mkdir ${SAMPLE}/${caller}/
        fi
        if [ ! -f "${SAMPLE}/${caller}/${SAMPLE}.${caller}.vcf" ]; then
            echo /global/home/users/pierrj/git/bash/run_${caller}.sh -f ${SAMPLE}/data/${SAMPLE}.bam -r ${REFERENCE} -t 1 -o ${SAMPLE}/${caller}/${SAMPLE}.${caller}.vcf -d ${SAMPLE}/${caller}/ >> ${caller}_jobqueue
        fi
        
    done < SAMPLE_LIST

    if [[ ! "${caller}" == "wham" ]]; then
        source activate /global/scratch/users/pierrj/conda-envs/${caller}
    fi

    parallel -j ${SLURM_NTASKS} < ${caller}_jobqueue

done