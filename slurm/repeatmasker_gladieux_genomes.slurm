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

cd /global/scratch/users/pierrj/moryzae_pav/te_annotations

RM_LIB=/global/scratch/users/pierrj/references/te_annotations/moryzae/repbase_plus_guy11_modeler_rm_database.fasta
GENOMES_PATH=gladieux_et_al_2021_assemblies
MAPFILE=gladieux_et_al_2021_assemblies_mapfile

module load repeatmasker/4.1.1s

while read genome; do
    RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${GENOMES_PATH}/${genome}_genomic.fna
done < ${MAPFILE}