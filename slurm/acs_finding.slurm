#!/bin/bash
#SBATCH --job-name=acs_finding
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/acs_finding

module load meme

if [ -d eccs_notes ]; then
    rm -r eccs_notes
fi

fimo -o eccs_notes acs_meme.txt notes.sorted.uniq.all.confirmedsplitreads.fasta