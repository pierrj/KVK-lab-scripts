#!/bin/bash
#SBATCH --job-name=space_sniffer_pav
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/

# du -h PAV_SV > pav_sv_du_02_28_24.txt

rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/full_run_tblastn_validation_2_7_2022/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/exonerate_validation/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/tree_making/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/tree_making_optimization/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/tree_problems/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/pav_validation/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/pav_validation_old_version/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/pav_validation/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast/pav_validation/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all/pav_validation/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast/orthofinder/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all/orthofinder/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/orthofinder/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/orthofinder_w_guy11/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/orthofinder_w_70_15/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/orthofinder_w_70_15/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/orthofinder_4_27_2022/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/tblastn_pav_validation_gladieux/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_anne_subset/orthofinder/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_anne_subset/pav_validation/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_anne_subset_re/orthofinder/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_anne_subset_re/pav_validation/
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/PAV_SV/PAV/adding_70_15_to_orthogrouping