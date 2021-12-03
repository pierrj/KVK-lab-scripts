#!/bin/bash
#SBATCH --job-name=repeatmasker_gladieux_genomes
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

# cd /global/scratch/users/pierrj/moryzae_pav/te_annotations

# RM_LIB=/global/scratch/users/pierrj/references/te_annotations/moryzae/repbase_plus_guy11_modeler_rm_database.fasta
# GENOMES_PATH=gladieux_et_al_2021_assemblies
# MAPFILE=gladieux_et_al_2021_assemblies_mapfile

# module load repeatmasker/4.1.1

# while read genome; do
#     RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${GENOMES_PATH}/${genome}.fasta
# done < ${MAPFILE}

# if [ -f "te_counts_per_genome" ]; then
#     rm te_counts_per_genome
# fi

# while read genome; do
#     line=$(wc -l ${genome}.fasta.out.gff | awk '{print $1}')
#     echo -e ${genome}'\t'${line} >> te_counts_per_genome
# done < ${MAPFILE}


cd /global/scratch/users/pierrj/moryzae_pav/te_vs_assembly_stats

GENOMES_PATH=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/
MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/mapfile

source activate /global/scratch/users/pierrj/conda_envs/gaas

export PERL5LIB=/global/scratch/users/pierrj/conda_envs/gaas/lib

while read genome; do
    gaas_fasta_statistics.pl -f ${GENOMES_PATH}/${genome}_genomic.fna | awk '{print $(NF-1)}' | awk '$1 !~ /-/' > ${genome}.stats
done < ${MAPFILE}