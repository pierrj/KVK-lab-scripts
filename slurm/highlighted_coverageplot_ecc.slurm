#!/bin/bash
#SBATCH --job-name=highlighted_coverageplot_splitreads
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
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina

genome_file=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes

while read sample
do
cd ${sample}
bedtools genomecov -d -i ${sample}.ecc_caller_out.splitreads.bed -g ${genome_file} > ${sample}.genomecoverage.bed
cd ..
done < mapfile

if [ -f "G3.normalization_table.total_splitreads" ]; then
    rm G3.normalization_table.total_splitreads
fi

while read sample; do
cd ${sample}
total=$(wc -l ${sample}.ecc_caller_out.splitreads.bed | awk '{print $1/100000}')
echo -e ${sample}'\t'${total} >> ../G3.normalization_table.total_splitreads
cd ..
done < mapfile

/global/home/users/pierrj/git/bash/generate_sample_biorep_treatment_mapfile_forme.sh -m mapfile

/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile -t G3_1A.genomecoverage.bed -n G3.normalization_table.total_splitreads -y t

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1 -b 100 -c 3 -n n

mapfile=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames
sample=G3
covfile=G3.normalized_binned
bedfile=/global/scratch/pierrj/references/te_annotations/moryzae/moryzae_copiaandgypsy_locs.bed

/global/home/users/pierrj/git/bash/generate_coverageplot_withhighlighting.sh -m ${mapfile} \
    -s ${sample} -c ${covfile} -b ${bedfile}

y_max=$(awk '{print $3}' ${sample}.normalized_binned.filtered | sort -nr | head -1 | awk '{print $0/1000}' | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}' | awk '{print $0*1000}')

echo ${y_max}