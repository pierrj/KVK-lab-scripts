#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
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
# for maggy here but copia1 is the same
LTR_FILE=maggy_sequence.fasta

if [ -f "ltr_circles_mapfile" ]; then
    rm ltr_circles_mapfile
fi

echo one_ltr_auto_circle_seq_randomized_maggy >> ltr_circles_mapfile
echo one_ltr_hr_circle_seq_randomized_maggy >> ltr_circles_mapfile
echo two_ltr_nhej_circle_seq_randomized_maggy >> ltr_circles_mapfile

bwa index ${LTR_FILE} -p ltr_db

while read circle; do
    # simulate paired end data 
    # 2k reads per scaffold, but 1000 scaffolds per file
    # for each circle sequence fasta
    /global/scratch/users/pierrj/art_bin_MountRainier/art_illumina -p -m 450 -s 50 -ss HSXt -i ${circle}.fasta -l 150 -c 2000 -o ${circle}.simulated
    echo "done simulating one"
    # map to the simple maggy sequence
    bwa mem -t ${THREADS} ltr_db ${circle}.simulated1.fq ${circle}.simulated2.fq -o ${circle}.sam
    samtools sort -O BAM ${circle}.sam > ${circle}.bam
    # calculate coverage across maggy sequence
    bedtools genomecov -ibam ${circle}.bam -d > ${circle}.coverage
done < ltr_circles_mapfile