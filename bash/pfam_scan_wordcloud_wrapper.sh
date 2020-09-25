#!/bin/bash
while getopts g:a:n:p:c:t: option
do
case "${option}"
in
g) SUBSET_GENE_IDS=${OPTARG};;
a) ALL_GENE_IDS=${OPTARG};;
n) OUTPUT_NAME=${OPTARG};;
p) PFAM_DIR=${OPTARG};;
c) CDS_FASTA=${OPTARG};;
t) THREADS=${OPTARG};;
esac
done

if [ -f "${OUTPUT_NAME}.subset.CDS.fasta" ]; then
    rm ${OUTPUT_NAME}.subset.CDS.fasta
fi

while read geneid;
do
grep -A1 gene=${geneid} ${CDS_FASTA} >> ${OUTPUT_NAME}.subset.CDS.fasta
done < ${SUBSET_GENE_IDS}

if [ -f "${OUTPUT_NAME}.all.CDS.fasta" ]; then
    rm ${OUTPUT_NAME}.all.CDS.fasta
fi

while read geneid;
do
grep -A1 gene=${geneid} ${CDS_FASTA} >> ${OUTPUT_NAME}.all.CDS.fasta
done < ${ALL_GENE_IDS}

pfam_scan.pl -cpu ${THREADS} -outfile ${OUTPUT_NAME}.subset.pfamscan.out -dir ${PFAM_DIR} -fasta ${OUTPUT_NAME}.subset.CDS.fasta

/global/scratch/users/pierrj/scripts/plant_rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl -p ${OUTPUT_NAME}.subset.pfamscan.out -e 0.001 -o ${OUTPUT_NAME}.subset.pfamscan.kparse.out

pfam_scan.pl -cpu ${THREADS} -outfile ${OUTPUT_NAME}.all.pfamscan.out -dir ${PFAM_DIR} -fasta ${OUTPUT_NAME}.all.CDS.fasta

/global/scratch/users/pierrj/scripts/plant_rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl -p ${OUTPUT_NAME}.all.pfamscan.out -e 0.001 -o ${OUTPUT_NAME}.all.pfamscan.kparse.out

python /global/home/users/pierrj/git/python/parse_pfam_scan_output.py ${OUTPUT_NAME}.subset.pfamscan.kparse.out ${OUTPUT_NAME}.all.pfamscan.kparse.out ${OUTPUT_NAME}.total.forwordcloud ${OUTPUT_NAME}.normalized.forwordcloud

Rscript --vanilla /global/home/users/pierrj/git/R/make_wordcloud.R ${OUTPUT_NAME}.total.forwordcloud ${OUTPUT_NAME}.total.wordcloud

Rscript --vanilla /global/home/users/pierrj/git/R/make_wordcloud.R ${OUTPUT_NAME}.normalized.forwordcloud ${OUTPUT_NAME}.normalized.wordcloud