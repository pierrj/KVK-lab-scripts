#!/bin/bash
#SBATCH --job-name=doesnt_matter
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd ${SLURM_SUBMIT_DIR}
cd ${sample}

/global/home/users/pierrj/git/bash/ena_generate_bam_file.sh -g ${genome_bwa} \
 -1 ${SLURM_SUBMIT_DIR}/${sample}/${sample}_1.fastq.gz -2 ${SLURM_SUBMIT_DIR}/${sample}/${sample}_2.fastq.gz \
 -s ${sample} -t ${SLURM_NTASKS} -m ${mapfile} -q ${READONE_LINK} -w ${READTWO_LINK}

/global/home/users/pierrj/git/bash/call_ecc_regions.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -b filtered.sorted.${sample}.bam

/global/home/users/pierrj/git/bash/assign_confidence_nodb.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -b filtered.sorted.${sample}.bam -r ${sample}.confirmedsplitreads.bed