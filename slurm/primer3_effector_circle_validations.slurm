#!/bin/bash
#SBATCH --job-name=primer3_effector_circle_validations
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


if [ -f "all_effectors_circles" ]; then
    rm all_effectors_circles
fi

while read effector_name effector_chrom effector_start effector_end; do
    if [ -f "${effector_name}_circles" ]; then
        rm ${effector_name}_circles
    fi
    while read sample; do
        awk -v OFS='\t' -v c=$effector_chrom -v s=$effector_start -v e=$effector_end '$1 == c && $2 < s && $3 > e' \
        /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.ecc_caller_out.details.nolowq.txt | \
            awk -v OFS='\t' -v s=$sample -v e=$effector_name '{print $1, $2, $3, e"_"s}' | awk -v OFS='\t' '{print $1, $2, $3, substr($4, length($4)-17, 17)}' >> ${effector_name}_circles
    done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
    ## pick one
    sort -k4,4 ${effector_name}_circles | head -1 >> all_effectors_circles
done < known_effector_locs


BEDFILE=all_effectors_circles
GENOME_FILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
OUTPUT_FILE=effector_circle_validation_primers.txt
NUM_PRIMER_PAIRS=5

/global/home/users/pierrj/git/bash/primer3_junction_primers.sh -b ${BEDFILE} -g ${GENOME_FILE} -o ${OUPUT_FILE} -n ${NUM_PRIMER_PAIRS}


MAKE SURE TO TRIPLE CHECK THAT PRIMERS LINE UP WITH THEIR PRIMER NAMES PLEASE