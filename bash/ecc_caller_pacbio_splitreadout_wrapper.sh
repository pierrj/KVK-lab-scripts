#!/bin/bash

## example wrapper script for ecc_caller_splitreadout.py

# get sam flags for reads
samtools view -F 4 ${effector_name}.${sample_target}.pacbioreads.bam | awk '{print $2}' > ${effector_name}.${sample_target}_aligned_samflags

# convert sam to bed and add back sam flags
bedtools bamtobed -cigar -i /global/scratch/users/pierrj/eccDNA/pipeline_tests/gviz_effector_ecc/${effector_name}.${sample_target}.pacbioreads.bam \
    | paste - ${effector_name}.${sample_target}_aligned_samflags > ${effector_name}.${sample_target}_aligned.bed

# call eccdnas and output split reads
python /global/home/users/pierrj/git/python/ecc_caller_pacbio_splitreadout.py ${effector_name}.${sample_target}_aligned.bed \
    ${effector_name}.${sample_target}_pacbio_eccs.bed 20 50000 ${effector_name}.${sample_target}_pacbio_splitreads.bed

# convert bed back to bam
~/.conda/envs/deeptools/bin/bedtools bedtobam -i ${effector_name}.${sample_target}_pacbio_splitreads.bed -g ${GENOME_CHROMSIZES} > ${effector_name}.${sample_target}_pacbio_splitreads.bam

samtools sort ${effector_name}.${sample_target}_pacbio_splitreads.bam > ${effector_name}.${sample_target}_pacbio_splitreads.sorted.bam

samtools index ${effector_name}.${sample_target}_pacbio_splitreads.sorted.bam