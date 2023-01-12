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
#MIT License
#
#Copyright (c) 2023 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

cd ${SLURM_SUBMIT_DIR}/${genome}

module purge

export PERL5LIB=''
export AUGUSTUS_CONFIG_PATH=/global/scratch/users/pierrj/fungap_runs/fungap_config_copy/
export AUGUSTUS_SCRIPTS_PATH=/opt/conda/bin

# run fungap through singularity
# using RNAseq data from SRA
# set up template directory for ease of use
singularity exec /global/scratch/users/pierrj/fungap_runs/fungap.sif.old python /workspace/FunGAP/fungap.py \
                                                        --output_dir fungap_out \
                                                        --trans_read_1 ../template_run/ERR5875670_1.fastq \
                                                        --trans_read_2 ../template_run/ERR5875670_2.fastq \
                                                        --genome_assembly ../genomes_to_annotate/${genome}_genomic.fna  \
                                                        --augustus_species magnaporthe_grisea  \
                                                        --sister_proteome ../template_run/prot_db.faa  \
                                                        --busco_dataset sordariomycetes_odb10 \
                                                        --num_cores ${SLURM_NTASKS}