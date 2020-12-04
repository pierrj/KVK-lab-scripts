#!/bin/bash
#SBATCH --job-name=eccDNA_gene_analysis_wrapper
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
source activate pfamscan

module remove r

module load r/3.4.2


cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/gene_analysis_wrapper/

while read treatment; do
    cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/gene_analysis_wrapper/${treatment}

    rm -r *

    cd /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress

    while read sample
    do
    cd ${sample}
    cp ecccaller_output.${sample}.renamed.details.tsv /global/scratch/users/pierrj/eccDNA/pipeline_tests/gene_analysis_wrapper/${treatment}
    cd ..
    done < mapfile_${treatment}

    cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/gene_analysis_wrapper/${treatment}

    PANNZER_INPUT=/global/scratch/users/pierrj/references/pannzer_GUY11_PacBio_merge_KROJ_post_repeatmasking.txt
    OUTPUT_NAME=${treatment}
    SCORE_CUTOFF=0.6
    MAPFILE=/global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/mapfile_${treatment}
    GENE_BEDFILE=/global/scratch/users/pierrj/references/GUY11_PacBio_merge_KROJ_renamed.just_genes.post_repeatmasking.T0name.bed
    TE_BEDFILE=/global/scratch/users/pierrj/references/GUY11_pacbio_ET_KROJ.features.justtranpsons.bed
    PFAM_DIR=/global/scratch/users/pierrj/pfam/
    CDS_FASTA=/global/scratch/users/pierrj/references/GUY11_PacBio_merge_KROJ_renamed.CDS.fasta
    THREADS=${SLURM_NTASKS}

    /global/home/users/pierrj/git/bash/gene_analysis_wrapper.sh -p ${PANNZER_INPUT} -n ${OUTPUT_NAME} -c ${SCORE_CUTOFF} -m ${MAPFILE} -b ${GENE_BEDFILE} -e ${TE_BEDFILE} -f ${PFAM_DIR} -a ${CDS_FASTA} -t ${THREADS}
done < treatment_mapfile