#!/bin/bash
#SBATCH --job-name=bs_test_only_some_scos
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/raxml_ng_test/only_some

# SCO_MSAS=/global/scratch/users/pierrj/PAV_SV/PAV/full_run_tree_1_31_2022/Single_Copy_Orthologue_Sequences_msas

# source activate /global/scratch/users/pierrj/conda_envs/orthofinder

# module purge

# ls -1 $SCO_MSAS > file_list

# shuf file_list | head -n 1000 > 1000_random_msas


# if [ -d "only_1000" ]; then
#     rm -r only_1000
# fi

# mkdir only_1000

# while read file; do
#     cp ${SCO_MSAS}/$file only_1000
# done < 1000_random_msas

# python /global/home/users/pierrj/git/python/concatenate_msas.py only_1000 ${ACCESSIONS} only_1000_msas_combined.fasta

# module load trimal

# trimal -gt 0.8 -in only_1000_msas_combined.fasta -out only_1000_msas_combined.trimmed.fasta

module purge
module load gcc/7.4.0
module load openmpi
module load cmake

# /global/scratch/users/pierrj/raxml_ng_savio1/bin/raxml-ng-mpi --parse --msa only_1000_msas_combined.trimmed.fasta --model PROTGTR+G --prefix only_1000

mpirun /global/scratch/users/pierrj/raxml_ng_savio1/bin/raxml-ng-mpi --bootstrap --bs-trees 1 --msa only_1000.raxml.rba --prefix only_1000 --seed 12345 --threads 10  --extra thread-pin