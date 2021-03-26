#!/bin/bash
#SBATCH --job-name=deletion_mummer_validation
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

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/deletion_mummer_plots/guy11_v_all/gene_deletion_output

source activate mummer
module load gnuplot
module load imagemagick

ls -1 | grep -v gene_mapfile > gene_mapfile

REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/guy11_genome_baoetal2017.fasta
GENOME_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes
SUBDIR_OUTPUT=output_jpgs
PERCENT_ZEROES_FILTER=1

while read gene; do
    cd ${gene}
    /global/home/users/pierrj/git/bash/auto_mummer_plot.sh -r ${REFERENCE} -g ${GENOME_DIR} -s ${SUBDIR_OUTPUT} -p ${PERCENT_ZEROES_FILTER}
    cd ..
done < gene_mapfile