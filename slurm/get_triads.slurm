#!/bin/bash
#SBATCH --job-name=get_triads
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/structural_variation_smoking_gun/anyoverlap

GFF_FILE=guy11_fungap_out_12_28_20.genes.gff3
ECC_MAPFILE=G3.eccmapfile_details

while read sample; do
    sample_base=$(basename ${sample})
    awk -v OFS='\t' '{if ($6 == "hconf") {print $1, $2, $3 }}' ${sample} > hconf.${sample_base}
    sort -k1,1 -k2,2n hconf.${sample_base} | awk -v OFS='\t' '{print $1, $2, $3, FNR}' > named_sorted.${sample_base}
    bedtools intersect -wa -c -a named_sorted.${sample_base} -b ${GFF_FILE} | awk '{if ($5 == 3){print $4}}' > triads_names.${sample_base}
    if [ -f "triads_splitreads.${sample_base}" ]; then
        rm triads_splitreads.${sample_base}
    fi
    while read splitread; do
        awk -v s=${splitread} '$4 == s' named_sorted.${sample_base} >> triads_splitreads.${sample_base}
    done < triads_names.${sample_base}
    bedtools intersect -F 1 -wa -wb -a triads_splitreads.${sample_base} -b ${GFF_FILE} | sort -k4,4n -k8,8n | awk '{print substr($13,4,10)".t1"}' > triads_genes.${sample_base}
    awk -v OFS='\t' '{if (NR%3==1) {
        first=$0
        getline
        second=$0
        getline
        {
            print first, second, $0
        }
    }
    }' triads_genes.${sample_base} | sort | uniq -c | awk -v OFS='\t' '{print $2, $3, $4, $1}' > triads_genes_counted.${sample_base}
done < ${ECC_MAPFILE}