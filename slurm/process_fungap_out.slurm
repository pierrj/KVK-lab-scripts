#!/bin/bash
#SBATCH --job-name=process_fungap_out
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

# cd /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast/

# if [ ! -d "all_proteomes_corrected" ]; then
#     mkdir all_proteomes_corrected
# fi

# source activate /global/scratch/users/pierrj/conda_envs/orthofinder
# module purge

# python /global/home/users/pierrj/git/python/process_protein_sequences_for_orthofinder.py all_proteomes \
#     all_proteomes_corrected \
#     gladieux_et_al_assembly_stats.txt \
#     output_accession_names

cd /global/scratch/users/pierrj/fungap_runs/wheat_blast/

MAPFILE=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_annotated_mapfile

while read genome; do
    grep -v '#' ${genome}/fungap_out/fungap_out/fungap_out.gff3 > ${genome}/fungap_out/fungap_out/fungap_out.nocomments.gff3
    sed -i '1s/^/##gff-version 3\n/' ${genome}/fungap_out/fungap_out/fungap_out.nocomments.gff3
done < $MAPFILE

/global/scratch/users/pierrj/conda_envs/orthofinder/bin/python /global/home/users/pierrj/git/python/process_gffs_for_orthofinder.py $MAPFILE

while read genome; do
    grep -v "#" all_gffs_fixed/${genome}_fungap_out.fixed.gff3 | grep -v gff-version=3=true | wc -l
done < $MAPFILE