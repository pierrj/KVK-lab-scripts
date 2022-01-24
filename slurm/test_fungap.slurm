#!/bin/bash
#SBATCH --job-name=test_fungap
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/fungap_runs/gladieux_all/test_run

module purge

export PERL5LIB=''
export AUGUSTUS_CONFIG_PATH=/global/scratch/users/pierrj/fungap_runs/fungap_config_copy/

## ADD GM_KEY TO HOME DIRECTORY ##

genome=ZJ2011-7-1

singularity exec --fakeroot --writable /global/scratch/users/pierrj/fungap_runs/fungap.sif python /workspace/FunGAP/fungap.py \
                                                        --output_dir fungap_out \
                                                        --trans_read_1 ERR5875670_1.fastq \
                                                        --trans_read_2 ERR5875670_2.fastq \
                                                        --genome_assembly ${genome}_genomic.fna  \
                                                        --augustus_species magnaporthe_grisea  \
                                                        --sister_proteome prot_db.faa  \
                                                        --busco_dataset sordariomycetes_odb10 \
                                                        --num_cores ${SLURM_NTASKS}