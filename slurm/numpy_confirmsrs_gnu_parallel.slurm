#!/bin/bash
#SBATCH --job-name=gnu_parallel
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/confirmsrs_numpy_gnu_parallel/
THREADS='20'
SAMPLE='CT_1A'
chrom_count='56'

split --number=l/${THREADS} --numeric-suffixes=1 --additional-suffix=.bed lengthfiltered.merged.splitreads.${SAMPLE}.renamed.bed lengthfiltered.merged.splitreads.${SAMPLE}.renamed.

parallel -j ${THREADS} --link python /global/home/users/pierrj/git/python/ecc_caller_anygenome_confirmsrs_numpy_gnuparallel.py lengthfiltered.merged.splitreads.${SAMPLE}.renamed.{}.bed sorted.grouped.outwardfacing.${SAMPLE}.renamed.bed ${SAMPLE} ${chrom_count} {} ::: $(seq -w 1 ${THREADS})

cat $(find . -maxdepth 1 -name "parallel.confirmed*" | xargs -r ls -1 | tr "\n" " ") > parallel.confirmed