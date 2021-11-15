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

# get effector circles

if [ -f "known_effector_locs" ]; then
    rm known_effector_locs
fi

# write location of effectors
echo -e AvrPita1'\t'MQOP01000011.1'\t'675417'\t'677485 >> known_effector_locs
echo -e AvrPita3'\t'MQOP01000016.1'\t'368546'\t'370280 >> known_effector_locs
echo -e AvrPi9'\t'MQOP01000008.1'\t'2277217'\t'2277733 >> known_effector_locs
echo -e AvrPi54'\t'MQOP01000009.1'\t'815984'\t'816445 >> known_effector_locs
echo -e AvrPiz-t'\t'MQOP01000012.1'\t'880675'\t'881001 >> known_effector_locs
echo -e Pwl4'\t'MQOP01000002.1'\t'4264949'\t'4265365 >> known_effector_locs

if [ -f "all_effectors_circles" ]; then
    rm all_effectors_circles
fi

# write all eccdnas containing effectors
while read effector_name effector_chrom effector_start effector_end; do
    if [ -f "${effector_name}_circles" ]; then
        rm ${effector_name}_circles
    fi
    while read sample; do
        # get eccdna locations first
        # then output with effector circle names
        awk -v OFS='\t' -v c=$effector_chrom -v s=$effector_start -v e=$effector_end '$1 == c && $2 < s && $3 > e' \
        /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.ecc_caller_out.details.nolowq.txt | \
            awk -v OFS='\t' -v s=$sample -v e=$effector_name '{print $1, $2, $3, e"_"s"ecc"}' | awk -v OFS='\t' '{print $1, $2, $3, substr($4, length($4)-17, 17)}' >> ${effector_name}_circles
    done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
    ## pick one effector circle to validate
    sort -k4,4 ${effector_name}_circles | head -1 >> all_effectors_circles
done < known_effector_locs

# exclude avrpita1 bc it only yields one good primer set which messes up primer3 out
grep -v 'AvrPita1' all_effectors_circles > all_effectors_nopita1_circles

# get fungicide resistance circles

if [ -f "fungicide_locs" ]; then
    rm fungicide_locs
fi

# fungicide gene locations
echo -e TRF1'\t'MQOP01000002.1'\t'2028651'\t'2031143 >> fungicide_locs
echo -e PTP2'\t'MQOP01000004.1'\t'3137138'\t'3138967 >> fungicide_locs

if [ -f "all_fungicide_circles" ]; then
    rm all_fungicide_circles
fi

# write all eccdnas containing fungicides
while read fungicide_name fungicide_chrom fungicide_start fungicide_end; do
    if [ -f "${fungicide_name}_circles" ]; then
        rm ${fungicide_name}_circles
    fi
    while read sample; do
        # get eccdna locations first
        # then output with effector circle names
        awk -v OFS='\t' -v c=$fungicide_chrom -v s=$fungicide_start -v e=$fungicide_end '$1 == c && $2 < s && $3 > e' \
        /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.ecc_caller_out.details.nolowq.txt | \
            awk -v OFS='\t' -v s=$sample -v e=$fungicide_name '{print $1, $2, $3, e"_"s"ecc"}' | awk -v OFS='\t' '{print $1, $2, $3, substr($4, length($4)-17, 17)}' >> ${fungicide_name}_circles
    done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
    ## pick one
    sort -k4,4 ${fungicide_name}_circles | head -1 >> all_fungicide_circles
done < fungicide_locs



BEDFILE=all_effectors_nopita1_circles
GENOME_FILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
OUTPUT_FILE=effector_circle_validation_nopita1_primers.txt
# primer pairs
NUM_PRIMER_PAIRS=2

## spit out primer pairs for to validate eccdna junctions
/global/home/users/pierrj/git/bash/primer3_junction_primers.sh -b ${BEDFILE} -g ${GENOME_FILE} -o ${OUTPUT_FILE} -n ${NUM_PRIMER_PAIRS}

BEDFILE=AvrPita1_circles
GENOME_FILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
OUTPUT_FILE=avrpita1_circle_validation_primers.txt
# primer pairs
NUM_PRIMER_PAIRS=1

## spit out primer pairs for to validate eccdna junctions
/global/home/users/pierrj/git/bash/primer3_junction_primers.sh -b ${BEDFILE} -g ${GENOME_FILE} -o ${OUTPUT_FILE} -n ${NUM_PRIMER_PAIRS}

BEDFILE=all_fungicide_circles
GENOME_FILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
OUTPUT_FILE=fungicide_circle_validation_primers.txt
# primer pairs
NUM_PRIMER_PAIRS=2

## spit out primer pairs for to validate eccdna junctions
/global/home/users/pierrj/git/bash/primer3_junction_primers.sh -b ${BEDFILE} -g ${GENOME_FILE} -o ${OUTPUT_FILE} -n ${NUM_PRIMER_PAIRS}

# bring it all together
cat effector_circle_validation_nopita1_primers.txt avrpita1_circle_validation_primers.txt fungicide_circle_validation_primers.txt > junction_validation_primers.txt
cat effector_circle_validation_nopita1_primers.txt.product_sizes avrpita1_circle_validation_primers.txt.product_sizes fungicide_circle_validation_primers.txt.product_sizes > junction_validation_primers.txt.product_sizes

cat all_effectors_nopita1_circles AvrPita1_circles all_fungicide_circles > all_validated_circle_locations.bed