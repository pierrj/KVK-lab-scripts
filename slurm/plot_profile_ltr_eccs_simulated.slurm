#!/bin/bash
#SBATCH --job-name=plot_profile_ltr_eccs_simulated
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/ltr_splitreads/simulated

THREADS=24
LTR_FILE=maggy_sequence.fasta

if [ -f "ltr_circles_mapfile" ]; then
    rm ltr_circles_mapfile
fi

echo one_ltr_auto_circle_seq_randomized_maggy >> ltr_circles_mapfile
echo one_ltr_hr_circle_seq_randomized_maggy >> ltr_circles_mapfile
echo two_ltr_nhej_circle_seq_randomized_maggy >> ltr_circles_mapfile

bwa index ${LTR_FILE} -p ltr_db

# probably 20k is the best but 2k temporarily
while read circle; do
    /global/scratch/users/pierrj/art_bin_MountRainier/art_illumina -p -m 450 -s 50 -ss HSXt -i ${circle}.fasta -l 150 -c 2000 -o ${circle}.simulated
    echo "done simulating one"
    bwa mem -t ${THREADS} ltr_db ${circle}.simulated1.fq ${circle}.simulated2.fq -o ${circle}.sam
    samtools sort -O BAM ${circle}.sam > ${circle}.bam
    bedtools genomecov -ibam ${circle}.bam -d > ${circle}.coverage
done < ltr_circles_mapfile