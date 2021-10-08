#!/bin/bash
#SBATCH --job-name=Pacbio_v_illumina_eccs
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pacbio

output_table=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pacbio/illumina_v_pacbio_out

if [ -f "${output_table}" ]; then
    rm ${output_table}
fi

while read sample; do
    python /global/home/users/pierrj/git/python/Pacbio_v_illumina_eccs.py /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pacbio/${sample}_pacbio_eccs.bed \
        /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.ecc_caller_out.details.nolowq.txt \
        /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pacbio/all_false_positives/${sample}/${sample}.mapq0_and_uniq_splitreads.allfalsepositives.bed \
        ${sample}.Pacbio_v_illumina ${sample}
        cat ${sample}.Pacbio_v_illumina >> ${output_table}
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile