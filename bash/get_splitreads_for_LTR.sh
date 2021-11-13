#!/bin/bash
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
while getopts s:b: option
do
case "${option}"
in
s) SAMPLE=${OPTARG};; # sample name
b) FILTERED_BAMFILE=${OPTARG};; # should be filtered to chromosomes/scaffolds of interest
esac
done

## USAGE ##
# this script filters reads from a bamfile to split reads that have potential to be attributed to ltr eccdnas
# see above for input descriptions

# get all split reads as long as they are more than 19 bp on either side and that they only appear exactly twice
samtools view -f 65 -F 4 ${FILTERED_BAMFILE} > tmp.filtered.sorted.allmapq.mapped.1.${SAMPLE}.sam
file='1'
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0}' tmp.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam

samtools view -f 129 -F 4 ${FILTERED_BAMFILE} > tmp.filtered.sorted.allmapq.mapped.2.${SAMPLE}.sam
file='2'
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0}' tmp.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam

# need separate files for read 1 and read 2, but put them together now
cat tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.1.${SAMPLE}.sam \
    tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.2.${SAMPLE}.sam > tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.all.${SAMPLE}.sam

# make sure that the split read length on one end matches the other
python /global/home/users/pierrj/git/python/filter_for_match_lengths.py tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.all.${SAMPLE}.sam tmp.match_length_filtered.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.all.${SAMPLE}.sam

# sort and index
samtools view -b -h <(cat <(samtools view -H ${FILTERED_BAMFILE}) tmp.match_length_filtered.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.all.${SAMPLE}.sam) > tmp.match_length_filtered.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.all.${SAMPLE}.bam
samtools sort tmp.match_length_filtered.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.all.${SAMPLE}.bam > candidate_ltr_srs.${SAMPLE}.bam
samtools index candidate_ltr_srs.${SAMPLE}.bam