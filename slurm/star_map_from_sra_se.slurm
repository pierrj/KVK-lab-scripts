#!/bin/bash
#SBATCH --job-name=doesntmatter
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

module load STAR

module load cufflinks

/global/home/users/pierrj/git/bash/star_map_from_sra_se.sh -g ${GENOME_DB} \
 -1 ${SLURM_SUBMIT_DIR}/${sample}/${sample}.sra.fastq \
 -s ${sample} -t ${SLURM_NTASKS}