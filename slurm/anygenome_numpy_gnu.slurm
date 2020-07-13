#!/bin/bash
#SBATCH --job-name=anygenome_numpy_gnu
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/anygenome_numpy_gnu

/global/home/users/pierrj/git/bash/call_ecc_regions.sh -m /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames\
    -s CT_1A -t 20 -b CT_1A.sorted.mergedandpe.bwamem.bam

/global/home/users/pierrj/git/bash/assign_confidence.sh -m /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames\
    -s CT_1A -t 20 -c CT_1A.genomecoverage.bed