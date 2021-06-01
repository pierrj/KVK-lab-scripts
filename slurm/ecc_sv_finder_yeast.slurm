#!/bin/bash
#SBATCH --job-name=ecc_sv_finder_yeast
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/structural_variation/yeast

module load gnuplot
module load imagemagick

cp /global/scratch/pierrj/eccDNA/pipeline_tests/test_ecc_sv_finder_script/yeast_test/EC1118_yeast.fasta .
cp /global/scratch/pierrj/eccDNA/pipeline_tests/test_ecc_sv_finder_script/yeast_test/M22_yeast.fasta .

/global/home/users/pierrj/git/bash/ecc_sv_finder.sh -r EC1118_yeast.fasta -q M22_yeast.fasta -t ${SLURM_NTASKS} -n ec1118_v_m22 \
    -p 1 -o ec1118_v_m22_out -m ec1118_v_m22_tmp


/global/home/users/pierrj/git/bash/ecc_sv_finder.sh -r M22_yeast.fasta -q EC1118_yeast.fasta -t ${SLURM_NTASKS} -n m22_v_ec1118 \
    -p 1 -o m22_v_ec1118_out -m m22_v_ec1118_tmp