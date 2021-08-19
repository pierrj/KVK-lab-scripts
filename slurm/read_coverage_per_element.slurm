#!/bin/bash
#SBATCH --job-name=read_coverage_per_element
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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/ltr_splitreads/reads_per_element_plot

TE_BEDFILE=/global/scratch/pierrj/references/te_annotations/moryzae/moryzae_high_conf_ltr_tes_nolengthfilter_individual_elements.bed
ELEMENT_MAPFILE=/global/scratch/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/coverage_plots/transposon_mapfile
MAPFILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
genome_file=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta

genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum;}' ${genome_file} )

while read ELEMENT; do
    grep ${ELEMENT} ${TE_BEDFILE} > ${ELEMENT}.bed
done < ${ELEMENT_MAPFILE}

if [ -f "te_eccs_per_sample.txt" ]; then
    rm te_eccs_per_sample.txt
fi

while read sample; do
    read_count=$(samtools view -c /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/no_secondary.filtered.sorted.${sample}.bam)
    while read ELEMENT; do
        bedtools coverage -sorted -a ${ELEMENT}.bed \
            -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/no_secondary.filtered.sorted.${sample}.bam | \
            awk -v OFS='\t' '{print $4, $5}' > ${ELEMENT}.${sample}.read_cov_perfeature
        element_bp=$(awk '{SUM+=$3-$2}END{print SUM}' ${ELEMENT}.bed)
        reads_per_element=$(awk '{SUM+=$2}END{print SUM}' ${ELEMENT}.${sample}.read_cov_perfeature)
        echo -e ${sample}'\t'${ELEMENT}'\t'${genome_size}'\t'${read_count}'\t'${element_bp}'\t'${reads_per_element} >> te_eccs_per_sample.txt
    done < ${ELEMENT_MAPFILE}
done < ${MAPFILE}