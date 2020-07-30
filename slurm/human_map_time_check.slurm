#!/bin/bash
#SBATCH --job-name=human_map_time_check
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/2018_moller/

mapfile="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.contignames"
genome_bwa="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GRCh37.p13_bwa"
sample="SRR6315407"
cd ${sample}

# date +"%T"

# /global/home/users/pierrj/git/bash/generate_bam_file_genomecov.sh -g ${genome_bwa} \
#  -1 ${SLURM_SUBMIT_DIR}/${sample}/${sample}_R1.fastq -2 ${SLURM_SUBMIT_DIR}/${sample}/${sample}_R2.fastq \
#  -s ${sample} -t ${SLURM_NTASKS}

# date +"%T"

# /global/home/users/pierrj/git/bash/generate_coverageplot.sh -m ${mapfile} -s ${sample} -c ${sample}.genomecoverage.bed

# date +"%T"

# /global/home/users/pierrj/git/bash/call_ecc_regions.sh -m ${mapfile} \
#     -s ${sample} -t ${SLURM_NTASKS} -b ${sample}.sorted.mergedandpe.bwamem.bam

date +"%T"

MAPFILE=${mapfile}
SAMPLE=${sample}
THREADS=${SLURM_NTASKS}
COVFILE=${sample}.genomecoverage.bed

# awk -v OFS='\t' 'NR==FNR{c[$1]++;next};c[$1]' ${MAPFILE} ${COVFILE} > ${SAMPLE}.genomecoverage.filtered.bed
# chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
# for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > tmp.chrom_count
# paste tmp.chrom_count ${MAPFILE} > tmp.chrom_count_and_names
# awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names ${SAMPLE}.genomecoverage.filtered.bed > ${SAMPLE}.genomecoverage.filtered.renamed.bed

# date +"%T"

# python /global/home/users/pierrj/git/python/merge_eccs.py ${SAMPLE} ${chrom_count}

# date +"%T"

# python /global/home/users/pierrj/git/python/make_coverage_db.py ${SAMPLE}.genomecoverage.filtered.renamed.bed ${chrom_count}

date +"%T"

split --number=l/${THREADS} --numeric-suffixes=1 merged.confirmed merged.confirmed

date +"%T"

# parallel -j ${THREADS} --link python /global/home/users/pierrj/git/python/coverage_confirm_db.py ${SAMPLE} {} ::: $(seq -w 1 ${THREADS})

python /global/home/users/pierrj/git/python/coverage_confirm_db.py ${SAMPLE} 

date +"%T"

cat $(find . -maxdepth 1 -name "ecccaller_output.${SAMPLE}.details.tsv*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.${SAMPLE}.details.tsv

cat $(find . -maxdepth 1 -name "ecccaller_output.${SAMPLE}.bed*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.${SAMPLE}.bed

paste ${MAPFILE} tmp.chrom_count > tmp.chrom_names_and_count

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.details.tsv > ecccaller_output.${SAMPLE}.renamed.bed

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.bed > ecccaller_output.${SAMPLE}.renamed.bed

date +"%T"