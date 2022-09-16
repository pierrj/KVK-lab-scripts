#!/bin/bash

genome=GCF_000001405.25_GRCh37.p13_genomic.fna
reads_1=SRR6315399.sra_1.fastq
reads_2=SRR6315399.sra_2.fastq

# params=default

# python ecc_finder.py map-sr ${genome} ${reads_1} ${reads_2} -r ${genome} -t 40 -o test_human >& run.test_human.${params}.out

# cp test_human/ecc.sr.csv test_params/${params}.csv

# rm test_human/ecc.sr.bam
# rm test_human/ecc.sr.bam.bed
# rm -r test_human/peak_files
# rm test_human/*



# params=peaks_100

# python ecc_finder.py map-sr ${genome} ${reads_1} ${reads_2} -r ${genome} -t 40 -d 100 -o test_human >& run.test_human.${params}.out

# cp test_human/ecc.sr.csv test_params/${params}.csv

# rm test_human/ecc.sr.bam
# rm test_human/ecc.sr.bam.bed
# rm -r test_human/peak_files
# rm test_human/*



# params=reads_0

# python ecc_finder.py map-sr ${genome} ${reads_1} ${reads_2} -r ${genome} -t 40 --min-read 0 -o test_human >& run.test_human.${params}.out

# cp test_human/ecc.sr.csv test_params/${params}.csv

# rm test_human/ecc.sr.bam
# rm test_human/ecc.sr.bam.bed
# rm -r test_human/peak_files
# rm test_human/*



# params=pvalue_1

# python ecc_finder.py map-sr ${genome} ${reads_1} ${reads_2} -r ${genome} -t 40 -p 1.0 -o test_human >& run.test_human.${params}.out

# cp test_human/ecc.sr.csv test_params/${params}.csv

# rm test_human/ecc.sr.bam
# rm test_human/ecc.sr.bam.bed
# rm -r test_human/peak_files
# rm test_human/*

params=peaks_10000

python ecc_finder.py map-sr ${genome} ${reads_1} ${reads_2} -r ${genome} -t 40 -d 10000 -o test_human >& run.test_human.${params}.out

cp test_human/ecc.sr.csv test_params/${params}.csv

rm test_human/ecc.sr.bam
rm test_human/ecc.sr.bam.bed
rm -r test_human/peak_files
rm test_human/*

params=length_0

python ecc_finder.py map-sr ${genome} ${reads_1} ${reads_2} -r ${genome} -t 40 -l 0 -o test_human >& run.test_human.${params}.out

cp test_human/ecc.sr.csv test_params/${params}.csv

rm test_human/ecc.sr.bam
rm test_human/ecc.sr.bam.bed
rm -r test_human/peak_files
rm test_human/*