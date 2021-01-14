#!/bin/bash
#SBATCH --job-name=get_multimapping_split_reads
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

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_calling_with_mapq0/retrieve_split_reads


SAMPLE=G3_1A
THREADS=20
MAPFILE="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
GENOME_DB="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"

date +"%T"

samtools sort -n ${SAMPLE}.sorted.mergedandpe.bwamem.bam > name_sorted.${SAMPLE}.sorted.mergedandpe.bwamem.bam

date +"%T"

samtools index name_sorted.${SAMPLE}.sorted.mergedandpe.bwamem.bam

date +"%T"

samtools view -b -F 256 name_sorted.${SAMPLE}.sorted.mergedandpe.bwamem.bam > primary_only.${SAMPLE}.sorted.mergedandpe.bwamem.bam

date +"%T"

/global/home/users/pierrj/git/bash/get_splitreads_for_LTR.sh -s ${SAMPLE} -b primary_only.${SAMPLE}.sorted.mergedandpe.bwamem.bam

date +"%T"

awk -v OFS='\t' '{
    prev=$0; f1=$1 ; f5=$5
    getline 
    if ($1 == f1 && $5 != 0 && f5 != 0) {
        print f1 > "tmp.doubleunique.readnames.1"
    }
    else if ($1 == f1 && $5 != 0 && f5 == 0) {
        print f1 > "tmp.singleunique.readnames.1"
    }
    else if ($1 == f1 && $5 == 0 && f5 != 0) {
        print f1 > "tmp.singleunique.readnames.1"
    }
    else if ($1 == f1 && $5 == 0 && f5 == 0) {
        print f1 > "tmp.doublemapq0.readnames.1"
    }
}'  tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.1.${SAMPLE}.sam

date +"%T"

awk -v OFS='\t' '{
    prev=$0; f1=$1 ; f5=$5
    getline 
    if ($1 == f1 && $5 != 0 && f5 != 0) {
        print f1 > "tmp.doubleunique.readnames.2"
    }
    else if ($1 == f1 && $5 != 0 && f5 == 0) {
        print f1 > "tmp.singleunique.readnames.2"
    }
    else if ($1 == f1 && $5 == 0 && f5 != 0) {
        print f1 > "tmp.singleunique.readnames.2"
    }
    else if ($1 == f1 && $5 == 0 && f5 == 0) {
        print f1 > "tmp.doublemapq0.readnames.2"
    }
}'  tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.2.${SAMPLE}.sam

date +"%T"

samtools view -b -f 65 -F 4 name_sorted.${SAMPLE}.sorted.mergedandpe.bwamem.bam > tmp.filtered.sorted.allmapq.mapped.1.${SAMPLE}.bam
samtools view -b -f 129 -F 4 name_sorted.${SAMPLE}.sorted.mergedandpe.bwamem.bam > tmp.filtered.sorted.allmapq.mapped.2.${SAMPLE}.bam

date +"%T"

java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar FilterSamReads \
        INPUT=tmp.filtered.sorted.allmapq.mapped.1.${SAMPLE}.bam \
        OUTPUT=${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bam \
        READ_LIST_FILE=tmp.singleunique.readnames.1 \
        FILTER=includeReadList \
        SORT_ORDER=unsorted

java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar FilterSamReads \
        INPUT=tmp.filtered.sorted.allmapq.mapped.2.${SAMPLE}.bam \
        OUTPUT=${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.2.bam \
        READ_LIST_FILE=tmp.singleunique.readnames.2 \
        FILTER=includeReadList \
        SORT_ORDER=unsorted

java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar FilterSamReads \
        INPUT=tmp.filtered.sorted.allmapq.mapped.1.${SAMPLE}.bam \
        OUTPUT=${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bam \
        READ_LIST_FILE=tmp.doublemapq0.readnames.1 \
        FILTER=includeReadList \
        SORT_ORDER=unsorted

java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar FilterSamReads \
        INPUT=tmp.filtered.sorted.allmapq.mapped.2.${SAMPLE}.bam \
        OUTPUT=${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.2.bam \
        READ_LIST_FILE=tmp.doublemapq0.readnames.2 \
        FILTER=includeReadList \
        SORT_ORDER=unsorted

date +"%T"

bedtools bamtobed -cigar -i ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bed
bedtools bamtobed -cigar -i ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.2.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.2.bed

date +"%T"

bedtools bamtobed -cigar -i ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bed
bedtools bamtobed -cigar -i ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.2.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.2.bed

date +"%T"

cat ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bed ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.2.bed > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.bed
cat ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bed ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.2.bed > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.bed

date +"%T"