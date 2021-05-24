#!/bin/bash
#!/bin/bash
while getopts m:s:t:b:q: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
b) FILTERED_BAMFILE=${OPTARG};; ## THIS SHOULD BE COORDINATE SORTED AND INCLUDE ONLY ONE RECORD PER ALIGNMENT AND IT SHOULD NOT INCLUDE MAPQ0 READS
q) FILTERED_BAMFILE_QSORTED=${OPTARG};; ## THIS SHOULD BE QNAME SORTED AND INCLUDE MULTIPLE RECORDS FOR MAPQ0 ALIGNMENTS
esac
done

samtools view -f 16 -F 4 ${FILTERED_BAMFILE} > tmp.reverseread1.${SAMPLE}.sam
splitread_file="reverseread1.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -F 20 ${FILTERED_BAMFILE} > tmp.forwardread1.${SAMPLE}.sam
splitread_file="forwardread1.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

cat tmp.oriented.samechromosome.exactlytwice.qualityfiltered.reverseread1.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.forwardread1.${SAMPLE}.sam > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam

python /global/home/users/pierrj/git/python/filter_for_match_lengths.py tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam tmp.match_length_filtered.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam

# converting to bed file
samtools view -b -h <(cat <(samtools view -H ${FILTERED_BAMFILE}) tmp.match_length_filtered.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam) > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.bam
bedtools bamtobed -i tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.bam | sort -k4,4 -k2,2n > splitreads.${SAMPLE}.bed

# merging split read halves into single, putative eccDNA forming regions to be confirmed or rejected
awk -v OFS='\t' '{
    prev=$0; f2=$2; f4=$4
    getline 
    if ($4 == f4 && f2 < $2) {
        print $1, f2, $3, $4
    }
}' splitreads.${SAMPLE}.bed > merged.splitreads.${SAMPLE}.bed

# length filter because we don't expect eccDNAs to be that big
# could be tweaked potentially but this gets rid of very few split reads
# 50k is based off the approximate size of eccDNAs that should be coming out of the column NEEDS TO BE VERIFIED
awk -v OFS='\t' '$3-$2<50000' merged.splitreads.${SAMPLE}.bed > lengthfiltered.merged.splitreads.${SAMPLE}.bed

# change names of scaffolds using mapfiles for compatability with any genome
chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > tmp.chrom_count
paste tmp.chrom_count ${MAPFILE} > tmp.chrom_count_and_names
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names lengthfiltered.merged.splitreads.${SAMPLE}.bed > lengthfiltered.merged.splitreads.${SAMPLE}.renamed.bed

awk -v OFS='\t' '{print $1, $2, $3}' lengthfiltered.merged.splitreads.${SAMPLE}.renamed.bed > unique_parallel.confirmed

awk '{print $3-$2}' unique_parallel.confirmed > dsn.unique_parallel.confirmed

samtools view -b -F 256 ${FILTERED_BAMFILE_QSORTED} > primary_only.${SAMPLE}.sorted.mergedandpe.bwamem.bam

samtools view -F 4 primary_only.${SAMPLE}.sorted.mergedandpe.bwamem.bam > tmp.primary_only.filtered.sorted.allmapq.mapped.1.${SAMPLE}.sam
file='1'
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0}' tmp.primary_only.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam

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

samtools view -b -F 4 ${FILTERED_BAMFILE_QSORTED} > tmp.filtered.sorted.allmapq.mapped.1.${SAMPLE}.bam

java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar FilterSamReads \
        INPUT=tmp.filtered.sorted.allmapq.mapped.1.${SAMPLE}.bam \
        OUTPUT=${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bam \
        READ_LIST_FILE=tmp.singleunique.readnames.1 \
        FILTER=includeReadList \
        SORT_ORDER=unsorted

java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar FilterSamReads \
        INPUT=tmp.filtered.sorted.allmapq.mapped.1.${SAMPLE}.bam \
        OUTPUT=${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bam \
        READ_LIST_FILE=tmp.doublemapq0.readnames.1 \
        FILTER=includeReadList \
        SORT_ORDER=unsorted

bedtools bamtobed -cigar -i ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bed
bedtools bamtobed -cigar -i ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bed

cp ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bed ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.bed
cp ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bed ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.bed

split --number=l/${THREADS} --numeric-suffixes=1 --additional-suffix=.bed ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.bed ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.chunk.

for i in $(seq -w 1 1 $THREADS); do
    python /global/home/users/pierrj/git/python/split_chunk_fixer.py ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.chunk.${i}.bed ${i}
done

seq -w 1 1 $((THREADS-1)) > tmp.seq
seq -w 2 1 ${THREADS} > tmp.seq_plusone
paste tmp.seq tmp.seq_plusone > tmp.seqs

while IFS=$'\t' read -r i next; do
    cat ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.chunk.${i}.bed split_line_fix.${next} > multimapped_splitreads.${i}.bed
done < tmp.seqs

cp ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.chunk.${THREADS}.bed multimapped_splitreads.${THREADS}.bed

parallel -j ${THREADS} --link python /global/home/users/pierrj/git/python/ecc_calling_mapq0.py dsn.unique_parallel.confirmed  multimapped_splitreads.{}.bed 50000 {} ::: $(seq -w 1 ${THREADS})

python /global/home/users/pierrj/git/python/ecc_calling_mapq0_singleunique.py dsn.unique_parallel.confirmed ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.bed 50000

cat mapq0_choices.* singleunique_choices > mapq0_single_unique_choices.bed

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names mapq0_single_unique_choices.bed > mapq0_single_unique_choices.renamed.bed

cp mapq0_single_unique_choices.renamed.bed mapq0_parallel.confirmed

cat unique_parallel.confirmed mapq0_parallel.confirmed > parallel.confirmed

paste ${MAPFILE} tmp.chrom_count > tmp.chrom_names_and_count
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count parallel.confirmed > ${SAMPLE}.confirmedsplitreads.bed

rm parallel.confirmed*
rm dsn.unique_parallel.confirmed
rm unique_parallel.confirmed
rm tmp.*
rm ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.*.bed
rm ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.bed
rm mapq0_choices.*
rm mapq0_single_unique_choices.bed
rm primary_only.${SAMPLE}.sorted.mergedandpe.bwamem.bam
rm lengthfiltered.merged.splitreads.${SAMPLE}.*
rm mapq0_parallel.confirmed
rm mapq0_single_unique_choices.renamed.*
rm multimapped_splitreads.*
rm merged.splitreads.${SAMPLE}.bed
rm split_line_fix.*