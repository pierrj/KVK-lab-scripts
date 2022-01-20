#!/bin/bash
#SBATCH --job-name=test_sv_calling_workshop_syri
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/sv_calling_workshop/syri_part

export CONDA_PKGS_DIRS=/global/scratch/users/pierrj/conda_envs/

source activate /global/scratch/users/pierrj/conda_envs/syri_2

# conda install -y python=3.8 cython=0.29.23 numpy=1.20.2 scipy=1.6.2 pandas=1.2.4 psutil=5.8.0 matplotlib=3.3.4
# conda install -y python-igraph=0.9.1
# conda install -y pysam=0.16.0.1
# conda install -y plotsr

# module purge

# git clone https://github.com/schneebergerlab/syri.git

# cd syri

# python setup.py install

# chmod +x syri/bin/syri syri/bin/chroder

module load mummer

nucmer --maxmatch -p c24_v_col-0_full -l 40 -g 90 -c 100 -b 200 -t ${SLURM_NTASKS} GCA_000001735.2_TAIR10.1_genomic.fna C24.chr.all.v2.0.fasta

delta-filter -m -i 90 -l 100 c24_v_col-0_full.delta > c24_v_col-0_full.filtered.delta

show-coords -THrd c24_v_col-0_full.filtered.delta > c24_v_col-0_full.filtered.coords

syri/syri/bin/syri -c c24_v_col-0_full.filtered.coords -d c24_v_col-0_full.filtered.delta -r GCA_000001735.2_TAIR10.1_genomic.fna -q C24.chr.all.v2.0.fasta --nc 5 --all -k

plotsr syri.out GCA_000001735.2_TAIR10.1_genomic.fna C24.chr.all.v2.0.fasta