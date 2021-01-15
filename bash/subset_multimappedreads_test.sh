
SAMPLE=G3_1A
THREADS=20
MAPFILE="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
GENOME_DB="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
FILTERED_BAMFILE=${SAMPLE}.sorted.mergedandpe.bwamem.bam

bwa mem -a -t ${THREADS} ${GENOME_DB} ${SAMPLE}_R1.sampled.fastq ${SAMPLE}_R2.sampled.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam

samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > ${SAMPLE}.mergedandpe.bwamem.bam
samtools sort -n ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

samtools view -b -F 256 ${SAMPLE}.sorted.mergedandpe.bwamem.bam > primary_only.${SAMPLE}.sorted.mergedandpe.bwamem.bam

/global/home/users/pierrj/git/bash/get_splitreads_for_LTR.sh -s ${SAMPLE} -b primary_only.${SAMPLE}.sorted.mergedandpe.bwamem.bam

## probably need proper bam files later

awk -v OFS='\t' '{
    prev=$0; f1=$1 ; f5=$5
    getline 
    if ($1 == f1 && $5 != 0 && f5 != 0) {
        print prev > "outfile1"
        print $0 > "outfile1"
    }
    else if ($1 == f1 && $5 != 0 && f5 == 0) {
        print prev > "outfile2"
        print $0 > "outfile2"
    }
    else if ($1 == f1 && $5 == 0 && f5 != 0) {
        print prev > "outfile2"
        print $0 > "outfile2"
    }
    else if ($1 == f1 && $5 == 0 && f5 == 0) {
        print prev > "outfile3"
        print $0 > "outfile3"
    }
}'  tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.1.${SAMPLE}.sam

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

samtools view -b -f 65 -F 4 ${FILTERED_BAMFILE} > tmp.filtered.sorted.allmapq.mapped.1.${SAMPLE}.bam
samtools view -b -f 129 -F 4 ${FILTERED_BAMFILE} > tmp.filtered.sorted.allmapq.mapped.2.${SAMPLE}.bam

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

bedtools bamtobed -cigar -i ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bed
bedtools bamtobed -cigar -i ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.2.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.2.bed

bedtools bamtobed -cigar -i ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bed
bedtools bamtobed -cigar -i ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.2.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.2.bed

cat ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bed ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.2.bed > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.bed
cat ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bed ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.2.bed > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.bed

## CHECK IF ALL SUPPLEMENTAL ALIGNMENTS HAVE MAPQ0



