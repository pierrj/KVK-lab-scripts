#!/bin/bash
#SBATCH --job-name=rm_proteoms
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/repeatmasked

RM_LIB=/global/scratch/users/pierrj/PAV_SV/PAV/te_annotations_anne/anne_repbase_te_lib_fixed
GENOMES_PATH=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/assemblies
MAPFILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/assemblies_mapfile

module load repeatmasker/4.1.1

echo test

while read genome; do
    RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${GENOMES_PATH}/${genome}
done < ${MAPFILE}