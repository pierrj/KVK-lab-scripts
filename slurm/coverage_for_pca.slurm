#!/bin/bash
#SBATCH --job-name=coverage_for_pca
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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pca_coverage_tech_reps_splitreads

CHROM_SIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes
WINDOWS=guy11_genome_baoetal2017.10000bpwindows
MAPFILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

bedtools makewindows -g ${CHROM_SIZES} -w 10000 > ${WINDOWS}

while read sample
do
    bedtools coverage -sorted -a ${WINDOWS} \
        -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.ecc_caller_out.splitreads.bed -g ${CHROM_SIZES} | \
        awk -v OFS='\t' '{print $4}' > ${sample}.genomecoverage
done < ${MAPFILE}


while read sample; do
    normalization=$(wc -l /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.ecc_caller_out.splitreads.bed | awk '{print $1/1000000}')
    echo -e ${sample}'\t'${normalization} >> normalization_table
done < ${MAPFILE}

while read sample normalization; do
    awk -v n=$normalization '{print $1/n}' ${sample}.genomecoverage > ${sample}.genomecoverage.normalized
done < normalization_table

awk -v OFS='_' '{print $1,$3}' ${WINDOWS} > first_column

paste first_column G3_1A.genomecoverage.normalized G3_1B.genomecoverage.normalized G3_1C.genomecoverage.normalized \
     G3_2A.genomecoverage.normalized G3_2B.genomecoverage.normalized G3_2C.genomecoverage.normalized \
      G3_3A.genomecoverage.normalized G3_3B.genomecoverage.normalized > all.genomecoverage.normalized

## transpose
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' all.genomecoverage.normalized > all.genomecoverage.normalized.transposed