#!/bin/bash
#SBATCH --job-name=acs_splitread_score
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

if [ -f split_read_scores ]; then
    rm split_read_scores
fi

while read scaffold start end; do
    awk -v c=$scaffold -v s=$start -v e=$end '{if ($2==c && $3==s && $4==e) print $1}' counted.uniq.all.confirmedsplitreads.bed >> split_read_scores
done < eccs_with_acs.new.bed