#!/bin/bash
#SBATCH --job-name=auto_mummer_plots_test
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/test_auto_plotting

source activate mummer
module load gnuplot
module load imagemagick

REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/guy11_genome_baoetal2017.fasta
REFERENCE_CHROMSIZE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/guy11_genome_baoetal2017.chromsizes
ISOLATE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/moryzae_70-15_ref.fasta
SUBDIR_INPUT=input_bed
SUBDIR_TMP=tmp_dir
SUBDIR_OUTPUT=output_jpgs

if [ -d "${SUBDIR_OUTPUT}" ]; then
    rm -r ${SUBDIR_OUTPUT}
fi

if [ -d "${SUBDIR_TMP}" ]; then
    rm -r ${SUBDIR_TMP}
fi

mkdir ${SUBDIR_OUTPUT}

mkdir ${SUBDIR_TMP}

alignment_number=$(ls -1 ${SUBDIR_INPUT} | wc -l | awk '{print ($1-1)/2}')

seq 0 ${alignment_number} > ${SUBDIR_TMP}/mapfile

while read plot; do
    bedtools getfasta -fi ${REFERENCE} -bed ${SUBDIR_INPUT}/guy11_alignment_${plot}.bed -fo ${SUBDIR_TMP}/guy11_alignment_${plot}.fasta
    samtools faidx ${SUBDIR_TMP}/guy11_alignment_${plot}.fasta
    cut -f1,2 ${SUBDIR_TMP}/guy11_alignment_${plot}.fasta.fai > ${SUBDIR_TMP}/guy11_alignment_${plot}.genomesize
    bedtools getfasta -fi ${ISOLATE} -bed ${SUBDIR_INPUT}/genome_alignment_${plot}.bed -fo ${SUBDIR_TMP}/genome_alignment_${plot}.fasta
    nucmer -p ${SUBDIR_TMP}/${plot} ${SUBDIR_TMP}/guy11_alignment_${plot}.fasta ${SUBDIR_TMP}/genome_alignment_${plot}.fasta
    show-coords ${SUBDIR_TMP}/${plot}.delta | tail -n +6 | awk -v OFS='\t' '{print $12, $1, $2}' | sort -k1,1 -k2,2n >  ${SUBDIR_TMP}/${plot}.bed
    bedtools genomecov -d -i ${SUBDIR_TMP}/${plot}.bed -g ${SUBDIR_TMP}/guy11_alignment_${plot}.genomesize > ${SUBDIR_TMP}/${plot}.genomecov
    total_size=$(wc -l ${SUBDIR_TMP}/${plot}.genomecov | awk '{print $1}')
    size_zeroes=$(awk '$3==0' ${SUBDIR_TMP}/${plot}.genomecov | wc -l | awk '{print $1}')
    percent_zeroes=$(awk -v var1=$size_zeroes -v var2=$total_size 'BEGIN { OFMT="%f";print  ( var1 / var2 ) }')
    if (( $(echo "$percent_zeroes < 0.1" |bc -l) ))
    then
        mummerplot --color -postscript -p ${SUBDIR_TMP}/${plot} ${SUBDIR_TMP}/${plot}.delta
        ps2pdf ${SUBDIR_TMP}/${plot}.ps ${SUBDIR_TMP}/${plot}.pdf
        convert -density 150 ${SUBDIR_TMP}/${plot}.pdf -quality 90 ${SUBDIR_OUTPUT}/${plot}.jpg
    fi
done < ${SUBDIR_TMP}/mapfile