#!/bin/bash
#SBATCH --job-name=run_fungap
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd ${SLURM_SUBMIT_DIR}/${genome}

module purge

export PERL5LIB=''
export AUGUSTUS_CONFIG_PATH=/global/scratch/users/pierrj/fungap_runs/fungap_config_copy/
export AUGUSTUS_SCRIPTS_PATH=/opt/conda/bin

singularity exec /global/scratch/users/pierrj/fungap_runs/fungap.sif.old python /workspace/FunGAP/fungap.py \
                                                        --output_dir fungap_out \
                                                        --trans_read_1 ../template_run/ERR5875670_1.fastq \
                                                        --trans_read_2 ../template_run/ERR5875670_2.fastq \
                                                        --genome_assembly ../Assembly/${genome}_final.scaffolds.fasta  \
                                                        --augustus_species magnaporthe_grisea  \
                                                        --sister_proteome ../template_run/prot_db.faa  \
                                                        --busco_dataset sordariomycetes_odb10 \
                                                        --num_cores ${SLURM_NTASKS}