#!/bin/bash
#SBATCH --job-name=validate_sanger
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pcr/

# process and label sanger sequencing 
cd sanger_11_8_21
    for i in *.seq; do
        sample_number=$(echo $i | awk '{print substr($1, 4, 4)}')
        echo $sample_number
        awk -v s=$sample_number 'BEGIN{print ">"s}; {print $0}' $i > $sample_number.fasta
    done
cd ..

if [ -f "sample_numbers" ]; then
    rm sample_numbers
fi

# write sample numbers
echo PJ06 >> sample_numbers
echo PJ09 >> sample_numbers
echo PJ18 >> sample_numbers
echo PJ24 >> sample_numbers
echo PJ30 >> sample_numbers
echo PJ33 >> sample_numbers
echo PJ39 >> sample_numbers
echo PJ45 >> sample_numbers

paste all_validated_circle_locations.bed sample_numbers > all_validated_circle_locations_w_sample_numbers.bed

GENOME_FILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta

## blast sanger sequence against fasta files of eccdnas
while read -r scaffold start end name sample_number; do
    # fasta file
    echo -e ${scaffold}'\t'${start}'\t'${end}'\t'${name} | bedtools getfasta -fi ${GENOME_FILE} -bed - -fo ${name}.fasta
    blastn -query sanger_11_8_21/${sample_number}.fasta -subject ${name}.fasta -out ${sample_number}_${name}.bast
done < all_validated_circle_locations_w_sample_numbers.bed

## print out blast results, some manual curation required
while read -r scaffold start end name sample_number; do
    echo ${sample_number}_${name}.bast
done < all_validated_circle_locations_w_sample_numbers.bed