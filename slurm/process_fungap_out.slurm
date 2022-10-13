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

# cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap

# if [ ! -d "all_proteomes" ]; then
#     mkdir all_proteomes
# fi

# MAPFILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/genomes_mapfile
# FUNGAP_DIR=/global/scratch/users/pierrj/fungap_runs/gladieux_all/

# while read genome; do
#     cp ${FUNGAP_DIR}/${genome}/fungap_out/fungap_out/fungap_out_prot.faa all_proteomes/${genome}_fungap_out_prot.faa
# done < $MAPFILE

# if [ ! -d "all_proteomes_corrected" ]; then
#     mkdir all_proteomes_corrected
# fi

# source activate /global/scratch/users/pierrj/conda_envs/orthofinder
# module purge

# /global/scratch/users/pierrj/conda_envs/orthofinder/bin/python /global/home/users/pierrj/git/python/process_protein_sequences_for_orthofinder.py all_proteomes \
#     all_proteomes_corrected \
#     gladieux_et_al_assembly_stats.txt \
#     output_accession_names

# /global/scratch/users/pierrj/conda_envs/orthofinder/bin/python /global/home/users/pierrj/git/python/process_protein_sequences_for_orthofinder.py all_proteomes_w_guy11 \
#     all_proteomes_corrected_w_guy11 \
#     gladieux_et_al_assembly_stats.txt \
#     output_accession_names_w_guy11

cd /global/scratch/users/pierrj/fungap_runs/wheat_blast

MAPFILE=/global/scratch/users/pierrj/fungap_runs/wheat_blast/wheat_blast_busco_greater_than_90_annotated

while read genome; do
    grep -v '#' ${genome}/fungap_out/fungap_out/fungap_out.gff3 > ${genome}/fungap_out/fungap_out/fungap_out.nocomments.gff3
    sed -i '1s/^/##gff-version 3\n/' ${genome}/fungap_out/fungap_out/fungap_out.nocomments.gff3
done < $MAPFILE

/global/scratch/users/pierrj/conda_envs/orthofinder/bin/python /global/home/users/pierrj/git/python/process_gffs_for_orthofinder.py $MAPFILE

while read genome; do
    mv all_gffs_fixed/${genome}_fungap_out.fixed.gff3 all_gffs_fixed/${genome}_fungap_out.fixed.gff3.old
    grep -v gff-version=3=true all_gffs_fixed/${genome}_fungap_out.fixed.gff3.old > all_gffs_fixed/${genome}_fungap_out.fixed.gff3
    rm all_gffs_fixed/${genome}_fungap_out.fixed.gff3.old
done < $MAPFILE