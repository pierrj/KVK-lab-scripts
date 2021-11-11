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
#SBATCH --job-name=profile_plots_readcov_per_element
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/ltr_splitreads/profileplot

# file of scaffold lengths
GENOME_CHROMSIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes
# mapfile of sample names
MAPFILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
# coverage file
READS_NAME=readcoverage
# reads file
READS_FILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/'${sample}'/no_secondary.filtered.sorted.'${sample}'.bam
WINDOWS=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/ltr_splitreads/profileplot/guy11_genome_baoetal2017.10bpwindows

source activate deeptools

bedtools makewindows -g ${GENOME_CHROMSIZES} -w 10 > guy11_genome_baoetal2017.10bpwindows

# eval to pass variables here
# get bedtools coverage file
while read sample
do
eval "bedtools coverage -sorted -a ${WINDOWS} \
    -b ${READS_FILE} -g ${GENOME_CHROMSIZES}" | \
    awk -v OFS='\t' '{print $1, $3, $4}' > ${sample}.${READS_NAME}
done < ${MAPFILE}

# generate sample mapfile
/global/home/users/pierrj/git/bash/generate_sample_biorep_treatment_mapfile_forme.sh -m ${MAPFILE}

# move samples
while read sample; do
    if [ -d "${sample}" ]; then
        rm -r ${sample}
    fi
    mkdir ${sample}
    mv ${sample}.${READS_NAME} ${sample}/${sample}.${READS_NAME}
done < ${MAPFILE}

if [ -f "${READS_NAME}.normalization_table" ]; then
    rm ${READS_NAME}.normalization_table
fi

# make normalization table
# based off number of reads in library
while read sample; do
    total=$( eval "samtools view -c ${READS_FILE}" | awk '{print $1/1000000}')
    echo -e ${sample}'\t'${total} >> ${READS_NAME}.normalization_table
done < ${MAPFILE}

# normalize and average read coverage across samples
/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile \
    -t G3_1A.${READS_NAME} -n ${READS_NAME}.normalization_table -y t

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1 -b 1 -c 3 -n n

mv G3.normalized_binned ${READS_NAME}.normalized_and_average

# convert to bedgraph
awk -v OFS='\t' '{print $1, $2-10, $2, $3}' ${READS_NAME}.normalized_and_average > ${READS_NAME}.normalized_and_average.bg

## get rid of overlapping entries and convert to bigwig
bedtools merge -i ${READS_NAME}.normalized_and_average.bg -c 4 -o mean -d -1 > ${READS_NAME}.normalized_and_average.corrected.bg
bedGraphToBigWig ${READS_NAME}.normalized_and_average.corrected.bg ${GENOME_CHROMSIZES} ${READS_NAME}.normalized_and_average.bw

ELEMENT_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/coverage_plots/transposon_mapfile
## high confidence, full length tes here from wickersoft, with length filter
TE_BEDFILE=/global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_high_conf_ltr_tes_for_plotprofile.bed

# generate profile pltos per element
while read ELEMENT; do
    grep ${ELEMENT} ${TE_BEDFILE} > ${ELEMENT}.bed
    if [[ $ELEMENT == GYMAG2 ]]
    then
        cp ${ELEMENT}.bed ${ELEMENT}.old.bed
        awk '$2 != 499211 && $3 != 505077' ${ELEMENT}.old.bed > ${ELEMENT}.bed ## remove one outlier entry for GYMAG2
    fi
    computeMatrix scale-regions -p ${SLURM_NTASKS} -S ${READS_NAME}.normalized_and_average.bw \
                            -R ${ELEMENT}.bed \
                            --beforeRegionStartLength 1000 \
                            --regionBodyLength 1000 \
                            --afterRegionStartLength 1000 \
                            --skipZeros -o ${ELEMENT}.${READS_NAME}.normalized_and_average.mat.gz ## calculate input for profile plots
    plotProfile -m ${ELEMENT}.${READS_NAME}.normalized_and_average.mat.gz \
            -out ${ELEMENT}.${READS_NAME}.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "${ELEMENT}" \
            --outFileNameData ${ELEMENT}.${READS_NAME}.normalized_and_average..tab ## output this to be plotted in ggplot
done < ${ELEMENT_MAPFILE}