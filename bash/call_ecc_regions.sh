#!/bin/bash
while getopts m:s:t:c:b: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
c) COVFILE=${OPTARG};;
b) SORTED_BAMFILE=${OPTARG};;
esac
done

samtools view -b -q 1 ${SORTED_BAMFILE} $(cat ${MAPFILE} | tr "\n" " ") > filtered.sorted.${SAMPLE}.bam

samtools view -f 81 -F 4 filtered.sorted.${SAMPLE}.bam > tmp.reverseread1.${SAMPLE}.sam
splitread_file="reverseread1.${SAMPLE}.sam"
awk '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if((a !~ /[DMIHS]/ && int(a) > 19 ) || (b !~ /[DMIHS]/ && int(b) > 19)) print $0}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -f 145 -F 4 filtered.sorted.${SAMPLE}.bam > tmp.reverseread2.${SAMPLE}.sam
splitread_file="reverseread2.${SAMPLE}.sam"
awk '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if((a !~ /[DMIHS]/ && int(a) > 19 ) || (b !~ /[DMIHS]/ && int(b) > 19)) print $0}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -f 65 -F 20 filtered.sorted.${SAMPLE}.bam > tmp.forwardread1.${SAMPLE}.sam
splitread_file="forwardread1.${SAMPLE}.sam"
awk '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if((a !~ /[DMIHS]/ && int(a) > 19 ) || (b !~ /[DMIHS]/ && int(b) > 19)) print $0}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -f 129 -F 20 filtered.sorted.${SAMPLE}.bam > tmp.forwardread2.${SAMPLE}.sam
splitread_file="forwardread2.${SAMPLE}.sam"
awk '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if((a !~ /[DMIHS]/ && int(a) > 19 ) || (b !~ /[DMIHS]/ && int(b) > 19)) print $0}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -f 16 -F 5 filtered.sorted.${SAMPLE}.bam > tmp.reversemerged.${SAMPLE}.sam
splitread_file="reversemerged.${SAMPLE}.sam"
awk '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if((a !~ /[DMIHS]/ && int(a) > 19 ) || (b !~ /[DMIHS]/ && int(b) > 19)) print $0}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -F 21 filtered.sorted.${SAMPLE}.bam > tmp.forwardmerged.${SAMPLE}.sam
splitread_file="forwardmerged.${SAMPLE}.sam"
awk '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if((a !~ /[DMIHS]/ && int(a) > 19 ) || (b !~ /[DMIHS]/ && int(b) > 19)) print $0}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

cat tmp.samechromosome.exactlytwice.qualityfiltered.reverseread1.${SAMPLE}.sam \
    tmp.samechromosome.exactlytwice.qualityfiltered.reverseread2.${SAMPLE}.sam \
    tmp.samechromosome.exactlytwice.qualityfiltered.forwardread1.${SAMPLE}.sam \
    tmp.samechromosome.exactlytwice.qualityfiltered.forwardread2.${SAMPLE}.sam \
    tmp.samechromosome.exactlytwice.qualityfiltered.reversemerged.${SAMPLE}.sam \
    tmp.samechromosome.exactlytwice.qualityfiltered.forwardmerged.${SAMPLE}.sam > tmp.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam

samtools view -b -h <(cat <(samtools view -H ${SORTED_BAMFILE}) tmp.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam) > tmp.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.bam
bedtools bamtobed -i tmp.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.bam | sort -k4,4 -k2,2n > splitreads.${SAMPLE}.bed

### NEED TO MERGE SPLIT READ ENTRIES HERE

### NEED TO FILTER SPLIT READ ENTRIES HERE

samtools view filtered.sorted.${SAMPLE}.bam | awk '{ if (($2 == 83 || $2 == 147 ) && $9 > 0) print $0 ; else if (($2 == 99 || $2 == 163) && $9 <0) print $0}' | cat <(samtools view -H ${SORTED_BAMFILE}) - | samtools view -b -h - > tmp.outwardfacing.${SAMPLE}.bam
bedtools bamtobed -i tmp.outwardfacing.${SAMPLE}.bam | sort -k 4,4 > tmp.outwardfacing.${SAMPLE}.bed
mv tmp.outwardfacing.${SAMPLE}.bed tmp.outwardfacing.${SAMPLE}.bed.old
awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,substr($4, 1, length($4)-2),$5,$6}' tmp.outwardfacing.${SAMPLE}.bed.old > tmp.outwardfacing.${SAMPLE}.bed.old.trimmed
awk 'NR==FNR{a[$4]++; next} a[$4]==2' tmp.outwardfacing.${SAMPLE}.bed.old.trimmed tmp.outwardfacing.${SAMPLE}.bed.old.trimmed > outwardfacing.${SAMPLE}.bed