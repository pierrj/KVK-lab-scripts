#!/bin/bash
while getopts m:s:b:t: option
do
case "${option}"
in
m) SPIKE_MAPFILE=${OPTARG};;
s) SAMPLE=${OPTARG};;
b) BAMFILE=${OPTARG};; ## needs index too
t) THREADS=${OPTARG};;
esac
done

if [ -f "spike_normalization_countable_${SAMPLE}" ]; then
    rm spike_normalization_countable_${SAMPLE}
fi
## input is a table with the ratio in column one and the spike in name as listed in the mapped fasta file in column two
while read ratio spike
do
echo ${spike} > tmp.mapfile_${spike}
samtools view -b -q 1 ${BAMFILE} $(cat tmp.mapfile_${spike} | tr "\n" " ") > tmp.filtered.sorted.${SAMPLE}_${spike}.bam
/global/home/users/pierrj/git/bash/call_ecc_regions.sh -m tmp.mapfile_${spike} \
    -s ${SAMPLE}_${spike} -t ${THREADS} -b tmp.filtered.sorted.${SAMPLE}_${spike}.bam
count=$(wc -l ${SAMPLE}_${spike}.confirmedsplitreads.bed | awk '{print $1}')
echo -e ${ratio}'\t'${spike}'\t'${count} >> spike_normalization_countable_${SAMPLE}
done < ${SPIKE_MAPFILE}

slope=$(awk '{ x[NR] = $1; y[NR] = $3;
 sx += x[NR]; sy += y[NR]; 
 sxx += x[NR]*x[NR];
 sxy += x[NR]*y[NR];
}
END{
 det = NR*sxx - sx*sx;
 a = (NR*sxy - sx*sy)/det;
 b = (-sx*sxy+sxx*sy)/det;
 print a/100000;
# for(i=1;i<=NR;i++) print x[i],a*x[i]+b;
}' spike_normalization_countable_${SAMPLE})

echo -e ${SAMPLE}'\t'${slope}