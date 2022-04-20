#!/bin/bash
#SBATCH --job-name=ecc_genome_cov_comparison
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/comparative/genome_coverage

if [ -f "ecc_coverage_comparison.txt" ]; then
    rm ecc_coverage_comparison.txt
fi

while read -r ORGANISM_NAME GENE_GFF SAMPLE_MAPFILE ECC_DIR GENOME_FILE COPIA_FILE GYPSY_FILE CONTIGNAMES_FILE; do
    # filter genome to contigs of interest, same as those used for ecc_caller
    samtools faidx ${GENOME_FILE} $(cat ${CONTIGNAMES_FILE} | tr "\n" " ") > ${GENOME_FILE}.filtered
    GENOME_FILE=$(realpath ${GENOME_FILE}.filtered)
    samtools faidx ${GENOME_FILE}
    cut -f1,2 ${GENOME_FILE}.fai > ${GENOME_FILE}.chromsizes
    GENOME_CHROMSIZES=$(realpath ${GENOME_FILE}.chromsizes)
    cd $ECC_DIR
    while read sample;
    do
    cd ${sample}
        awk ' $5 != "lowq" ' ${sample}.ecc_caller_out.details.txt > ${sample}.ecc_caller_out.details.nolowq.txt
        ## base pairs
        genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum;}' ${GENOME_FILE} )
        ## coverage base pairs
        covered_bps=$(bedtools genomecov -d -i ${sample}.ecc_caller_out.details.nolowq.txt -g ${GENOME_FILE}.chromsizes | awk '$3 == 0' | wc -l)
        echo -e $sample'\t'$genome_size'\t'$covered_bps >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/comparative/genome_coverage/ecc_coverage_comparison.txt
    cd ..
    done < $SAMPLE_MAPFILE
done < mapfile_ecc_count

paste ecc_coverage_comparison.txt ecc_count_metadata.txt > ecc_coverage_comparison_w_metadata.txt