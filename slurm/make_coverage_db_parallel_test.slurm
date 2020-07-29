#!/bin/bash
#SBATCH --job-name=make_coverage_db_parallel_test
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/makecoverage_db_parallel

MAPFILE="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
SAMPLE=CT_1A
THREADS=${SLURM_NTASKS}
COVFILE=${sample}.genomecoverage.bed

awk -v OFS='\t' 'NR==FNR{c[$1]++;next};c[$1]' ${MAPFILE} ${COVFILE} > ${SAMPLE}.genomecoverage.filtered.bed
chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > tmp.chrom_count
paste tmp.chrom_count ${MAPFILE} > tmp.chrom_count_and_names
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names ${SAMPLE}.genomecoverage.filtered.bed > ${SAMPLE}.genomecoverage.filtered.renamed.bed

date +"%T"

for (( i = 1 ; i < ${chrom_count}+1; i++))
do 
    awk -v I=$i '$1==I' ${SAMPLE}.genomecoverage.filtered.renamed.bed > ${SAMPLE}.$i.genomecoverage.filtered.renamed.bed
done

parallel -j ${chrom_count} --link python /global/home/users/pierrj/git/python/make_coverage_db_parallel.py ${SAMPLE}.{}.genomecoverage.filtered.renamed.bed {} ::: $(seq ${chrom_count})

date +"%T"

split --number=l/${THREADS} --numeric-suffixes=1 merged.confirmed merged.confirmed
parallel -j ${THREADS} --link python /global/home/users/pierrj/git/python/coverage_confirm_db.py ${SAMPLE} {} ::: $(seq -w 1 ${THREADS})
cat $(find . -maxdepth 1 -name "ecccaller_output.${SAMPLE}.details.tsv*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.${SAMPLE}.details.tsv

cat $(find . -maxdepth 1 -name "ecccaller_output.${SAMPLE}.bed*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.${SAMPLE}.bed

paste ${MAPFILE} tmp.chrom_count > tmp.chrom_names_and_count

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.details.tsv > ecccaller_output.${SAMPLE}.renamed.bed

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.bed > ecccaller_output.${SAMPLE}.renamed.bed

md5sum ecccaller_output.${SAMPLE}.renamed.bed
md5sum ecccaller_output.${SAMPLE}.renamed.bed

ls | grep -v CT_1A.genomecoverage.bed | xargs rm

date +"%T"

python /global/home/users/pierrj/git/python/make_coverage_db.py ${SAMPLE}.genomecoverage.filtered.renamed.bed ${chrom_count}

date +"%T"

split --number=l/${THREADS} --numeric-suffixes=1 merged.confirmed merged.confirmed
parallel -j ${THREADS} --link python /global/home/users/pierrj/git/python/coverage_confirm_db.py ${SAMPLE} {} ::: $(seq -w 1 ${THREADS})
cat $(find . -maxdepth 1 -name "ecccaller_output.${SAMPLE}.details.tsv*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.${SAMPLE}.details.tsv

cat $(find . -maxdepth 1 -name "ecccaller_output.${SAMPLE}.bed*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.${SAMPLE}.bed

paste ${MAPFILE} tmp.chrom_count > tmp.chrom_names_and_count

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.details.tsv > ecccaller_output.${SAMPLE}.renamed.bed

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.bed > ecccaller_output.${SAMPLE}.renamed.bed

md5sum ecccaller_output.${SAMPLE}.renamed.bed
md5sum ecccaller_output.${SAMPLE}.renamed.bed

ls | grep -v CT_1A.genomecoverage.bed | xargs rm

date +"%T"