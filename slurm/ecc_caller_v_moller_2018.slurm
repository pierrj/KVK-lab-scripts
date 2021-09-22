#!/bin/bash
#SBATCH --job-name=ecc_caller_v_moller_2018
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/human

sample_mapfile=/global/scratch/users/pierrj/eccDNA/2018_moller/full_run/mapfile
mapfile="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.contignames"
muscle_mapfile=/global/scratch/users/pierrj/eccDNA/2018_moller/full_run/mapfile_muscle_moller_samples
leukocytes_mapfile=/global/scratch/users/pierrj/eccDNA/2018_moller/full_run/mapfile_leukocytes_moller_samples
moller_muscle_eccs=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/human/2018_moller_muscle_eccs.txt
moller_leukocyte_eccs=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/human/2018_moller_leukocyte_eccs.txt
output_table=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/human/ecc_caller_v_moller_2018_out_table.txt

# while read sample; do
#     if [ -d "${sample}" ]; then
#         rm -r ${sample}
#     fi
#     mkdir ${sample}
#     cd ${sample}
#         # get uniquely mapped eccdnas
#         /global/home/users/pierrj/git/bash/ecc_caller_only_unique.sh -m ${mapfile} \
#             -s ${sample} \
#             -t ${SLURM_NTASKS} \
#             -b /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${sample}/uniq.filtered.sorted.${sample}.bam
#         /global/home/users/pierrj/git/bash/assign_confidence_nodb_nomerge.sh -m ${mapfile} \
#             -s ${sample} \
#             -t ${SLURM_NTASKS} \
#             -b /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${sample}/no_secondary.filtered.sorted.${sample}.bam \
#             -r ${sample}.unique.confirmedsplitreads.bed
#         cp ${sample}.ecc_caller_out.details.txt ${sample}.ecc_caller_out.uniq.details.txt
#         # get all eccdnas
#         cp /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/${sample}/${sample}.ecc_caller_out.details.txt ${sample}.ecc_caller_out.details.txt
#     cd ..
# done < ${sample_mapfile}

if [ -f "${output_table}" ]; then
    rm ${output_table}
fi

while read sample moller_sample; do
    cd ${sample}
    python /global/home/users/pierrj/git/python/ecc_caller_v_moller_2018.py ${moller_muscle_eccs} ${moller_sample} ${sample}.ecc_caller_out.details.txt ${sample}.ecc_caller_out.uniq.details.txt muscle ${sample} ${sample}.ecc_caller_v_moller_2018
    cat ${sample}.ecc_caller_v_moller_2018 >> ${output_table}
    cd ..
done < ${muscle_mapfile}



while read sample moller_sample; do
    cd ${sample}
    echo stuff?
    python /global/home/users/pierrj/git/python/ecc_caller_v_moller_2018.py ${moller_leukocyte_eccs} ${moller_sample} ${sample}.ecc_caller_out.details.txt ${sample}.ecc_caller_out.uniq.details.txt leukocytes ${sample} ${sample}.ecc_caller_v_moller_2018
    cat ${sample}.ecc_caller_v_moller_2018 >> ${output_table}
    cd ..
done < ${leukocytes_mapfile}