#!/bin/bash
#SBATCH --job-name=auto_mummer_plots_test
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/test_auto_plotting

conda activate mummer

REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/guy11_genome_baoetal2017.fasta
ISOLATE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/moryzae_70-15_ref.fasta
SUBDIR=test_plot_script
SUBDIR_OUTPUT=test_plot_output

while read plot; do
    bedtools getfasta -fi ${REFERENCE} -bed ${SUBDIR}/guy11_alignment_${plot}.bed -fo ${SUBDIR}/guy11_alignment_${plot}.fasta
    bedtools getfasta -fi ${ISOLATE} -bed ${SUBDIR}/genome_alignment_${plot}.bed -fo ${SUBDIR}/genome_alignment_${plot}.fasta
    nucmer -p ${SUBDIR}/${plot} ${SUBDIR}/guy11_alignment_${plot}.fasta ${SUBDIR}/genome_alignment_${plot}.fasta
    mummerplot --color -postscript -p ${SUBDIR}/${plot} ${SUBDIR}/${plot}.delta
    ps2pdf ${SUBDIR}/${plot}.ps ${SUBDIR_OUTPUT}/${plot}.pdf
done < ${SUBDIR}/mapfile