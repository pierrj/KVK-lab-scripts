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

samtools view -F 4 ${SORTED_BAMFILE} $(cat ${MAPFILE} | tr "\n" " ") > tmp.filtered.sorted.allmapq.mapped.${SAMPLE}.sam
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0}' tmp.filtered.sorted.allmapq.mapped.${SAMPLE}.sam > tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${SAMPLE}.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${SAMPLE}.bam tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${SAMPLE}.sam > tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.${SAMPLE}.sam

samtools view -b -h <(cat <(samtools view -H ${SORTED_BAMFILE}) tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.${SAMPLE}.sam) > tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.${SAMPLE}.bam

samtools sort tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.${SAMPLE}.bam > SRs.no_orientation.allmapq.bam
samtools index SRs.no_orientation.allmapq.bam