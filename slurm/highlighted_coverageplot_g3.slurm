#!/bin/bash
#SBATCH --job-name=highlighted_coverageplot_g3
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

while read sample
do
cd ${sample}
bedtools genomecov -d -ibam ${sample}.sorted.mergedandpe.bwamem.bam > ${sample}.genomecoverage.bed
cd ..
done < mapfile

/global/home/users/pierrj/git/bash/generate_sample_biorep_treatment_mapfile_forme.sh -m mapfile

/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile -t G3_1A.genomecoverage.bed -n filtered.sorted.G3_1A.bam -y e

mv mapfile_for_normalize_and_average mapfile_for_normalize_and_average_guy11reads

cd results/coverageplot

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ../../mapfile_for_normalize_and_average_guy11reads -f 1000000 -b 100 -c 3 -n f

mapfile=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames
sample=G3
covfile=G3.normalized_binned
bedfile=/global/scratch/users/pierrj/repeatmasker/moryzae_ltr_parts/guy11_genome_baoetal2017.ltr_parts.bed

/global/home/users/pierrj/git/bash/generate_coverageplot_withhighlighting.sh -m ${mapfile} \
    -s ${sample} -c ${covfile} -b ${bedfile}

y_max=$(awk '{print $3}' ${SAMPLE}.normalized_binned.filtered | sort -nr | head -1 | awk '{print $0/1000}' | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}' | awk '{print $0*1000}')

echo ${y_max}