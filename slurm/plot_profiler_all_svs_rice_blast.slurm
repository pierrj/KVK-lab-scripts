#!/bin/bash
#SBATCH --job-name=plot_profiler_all_svs_rice_blast
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
#Copyright (c) 2023 Pierre Michel Joubert
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

WORKING_DIR=/global/scratch/users/pierrj/PAV_SV/SV/dels_re_rice_blast

cd $WORKING_DIR

FASTA_FILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta

basename_fasta_file=$(basename $FASTA_FILE)

SV_DIR=/global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11
SV_MAPFILE=gladieux_accessions

## get svs
if [ -f all_svs_list ]; then
    rm all_svs_list
fi
while read SAMPLE; do
    echo ${SV_DIR}/${SAMPLE}/all/${SAMPLE}.all.vcf >> all_svs_list
done < ${SV_DIR}/${SV_MAPFILE}

## survivor merge positional args
# in order #
# distance between breakpoints
# number of supporting callers
# take type into account
# take strand into account
# disabled?
# minimum size
# output
/global/scratch/users/pierrj/SURVIVOR/Debug/SURVIVOR merge all_svs_list 100 1 1 1 0 0 all_svs.vcf
/global/scratch/users/pierrj/SURVIVOR/Debug/SURVIVOR vcftobed all_svs.vcf 1 10000000000000000 all_svs.bed

## remove mitochondria from svs
mv all_svs.bed all_svs.bed.old
grep -v 'Supercontig_7.9' all_svs.bed.old > all_svs.bed


## get TE annotation file
if [ ! -f "${basename_fasta_file}.out.gff" ]; then
    module load repeatmasker/4.1.1
    RM_LIB=/global/scratch/users/pierrj/PAV_SV/PAV/te_annotations_anne/anne_repbase_te_lib_fixed
    RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . $FASTA_FILE
fi

source activate /global/scratch/users/pierrj/conda_envs/deeptools

## generate some necessary files
samtools faidx ${FASTA_FILE}
cut -f1-2 ${FASTA_FILE}.fai > ${FASTA_FILE}.chromsizes

sv=DEL
GENOME_FILE=$FASTA_FILE
GENOME_CHROMSIZES=${FASTA_FILE}.chromsizes
THREADS=${SLURM_NTASKS}
WINDOWS=10
GENE_BED=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11/GUY11.genes.bed
TE_GFF=${basename_fasta_file}.out.gff

# get subset files for tes, gene deletions
awk -v sv=$sv -v OFS='\t' '{ if ( $11 == sv ) {print $12, $13, $15}}' all_svs.bed > all_${sv}.bed
bedtools intersect -a all_${sv}.bed -b $GENE_BED -wa -u -f 0.5 > gene_${sv}.bed
bedtools intersect -a all_${sv}.bed -b $TE_GFF -wa -u -f 0.5 > te_${sv}.bed

# replicate many times to get genomic baseline
for i in $(seq 1 100); do cat all_${sv}.bed; done > all_${sv}.bed.replicated

bedtools shuffle -i all_${sv}.bed.replicated -g $GENOME_CHROMSIZES > shuffled_${sv}.bed

awk -v OFS='\t' '{print $1, $4, $5}' $TE_GFF | grep -v "#" > all_tes.bed

cp $GENE_BED all_genes.bed

## calculate te and gene density for all deletion subsets and genomic baseline
for REGIONS_BED in all_${sv}.bed gene_${sv}.bed te_${sv}.bed shuffled_${sv}.bed all_genes.bed all_tes.bed; do
    OUTPUT_NAME=plot_profile_${REGIONS_BED}

    DENSITY_BED=$TE_GFF
    SUB_OUTPUT_NAME=tes
    FLANKING_DIST=5000
    /global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv} -f $FLANKING_DIST

    DENSITY_BED=$GENE_BED
    SUB_OUTPUT_NAME=genes
    FLANKING_DIST=5000
    /global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv} -f $FLANKING_DIST
done