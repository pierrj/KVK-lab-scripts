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
cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pcr/

# cd sanger_11_8_21
#     for i in *.seq; do
#         sample_number=$(echo $i | awk '{print substr($1, 4, 4)}')
#         echo $sample_number
#         awk -v s=$sample_number 'BEGIN{print ">"s}; {print $0}' $i > $sample_number.fasta
#     done
# cd ..

if [ -f "sample_numbers" ]; then
    rm sample_numbers
fi

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

while read -r scaffold start end name sample_number; do
    echo -e ${scaffold}'\t'${start}'\t'${end}'\t'${name} | bedtools getfasta -fi ${GENOME_FILE} -bed - -fo ${name}.fasta
    blastn -query sanger_11_8_21/${sample_number}.fasta -subject ${name}.fasta -out ${sample_number}_${name}.bast
done < all_validated_circle_locations_w_sample_numbers.bed


while read -r scaffold start end name sample_number; do
    echo ${sample_number}_${name}.bast
done < all_validated_circle_locations_w_sample_numbers.bed