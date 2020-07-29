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

MAPFILE="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.contignames"
SAMPLE=SRR6315407
THREADS=${SLURM_NTASKS}

chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')

date +"%T"

for (( i = 1 ; i < ${chrom_count}+1; i++))
do 
    awk -v I=$i '$1==I' ${SAMPLE}.genomecoverage.filtered.renamed.bed > ${SAMPLE}.$i.genomecoverage.filtered.renamed.bed
done

parallel -j ${chrom_count} --link python /global/home/users/pierrj/git/python/make_coverage_db_parallel.py ${SAMPLE}.{}.genomecoverage.filtered.renamed.bed {} ::: $(seq ${chrom_count})

date +"%T"

ls | grep -v ${SAMPLE}.genomecoverage.filtered.renamed.bed | xargs rm

date +"%T"

python /global/home/users/pierrj/git/python/make_coverage_db.py ${SAMPLE}.genomecoverage.filtered.renamed.bed ${chrom_count}

date +"%T"

ls | grep -v ${SAMPLE}.genomecoverage.filtered.renamed.bed | xargs rm

date +"%T"