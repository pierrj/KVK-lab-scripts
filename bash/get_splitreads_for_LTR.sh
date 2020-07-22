#!/bin/bash
while getopts m:s:t:b: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
s) SAMPLE=${OPTARG};;
b) SORTED_BAMFILE=${OPTARG};; ### THE INDEX FILE HAS TO BE HERE TOO
esac
done

samtools view -b ${SORTED_BAMFILE} $(cat ${MAPFILE} | tr "\n" " ") > tmp.filtered.sorted.allmapq.${SAMPLE}.bam

samtools view -f 65 -F 4 tmp.filtered.sorted.allmapq.${SAMPLE}.bam > tmp.filtered.sorted.allmapq.mapped.1.${SAMPLE}.sam
file='1'
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0}' tmp.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam

samtools view -f 129 -F 4 tmp.filtered.sorted.allmapq.${SAMPLE}.bam > tmp.filtered.sorted.allmapq.mapped.2.${SAMPLE}.sam
file='2'
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0}' tmp.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam

samtools view -F 5 tmp.filtered.sorted.allmapq.${SAMPLE}.bam > tmp.filtered.sorted.allmapq.mapped.merged.${SAMPLE}.sam
file='merged'
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0}' tmp.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam

cat tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.1.${SAMPLE}.sam \
    tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.2.${SAMPLE}.sam \
    tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.merged.${SAMPLE}.sam > tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.all.${SAMPLE}.sam

samtools view -b -h <(cat <(samtools view -H ${SORTED_BAMFILE}) tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.all.${SAMPLE}.sam) > tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.all.${SAMPLE}.bam

samtools sort tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.all.${SAMPLE}.bam > SRs.no_orientation.allmapq.${SAMPLE}.bam
samtools index SRs.no_orientation.allmapq.${SAMPLE}.bam