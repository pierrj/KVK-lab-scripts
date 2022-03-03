#!/bin/bash
#SBATCH --job-name=samtools_stats
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

# cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio/SMRTcell_1

# samtools stats SMRTcell_1_subreads.bam
# # 23329067522

# samtools stats SMRTcell_1_ccs.bam
# # 1242896329

# cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio/SMRTcell_2

# samtools stats SMRTcell_2_subreads.bam
# # 25008001988
# samtools stats SMRTcell_2_ccs.bam
# # 1365935152

# cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio/SMRTcell_3

# samtools stats SMRTcell_3_subreads.bam
# # 23170453445

# samtools stats SMRTcell_3_ccs.bam
# # 1307942209

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read ${sample}; do
    echo ${sample}
    cat ${sample}/${sample}_R1.fastq | paste - - - - | cut -f 2 | tr -d '\n' | wc -c
    cat ${sample}/${sample}_R2.fastq | paste - - - - | cut -f 2 | tr -d '\n' | wc -c
done < mapfile