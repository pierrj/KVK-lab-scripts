
sample=ERR2660591
sample=SRR16282278
sample=SRR11528297


awk -v c=$effector_chrom -v s=$region_start -v e=$region_end '$1 == c && $2 > s && $3 < e' ${sample}/${sample}.ecc_caller_out.details.txt


GENE_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.renamed.bed
ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/neverfound/G3.all_eccdnas_mapfile

5,6,7

ECCDNA_FILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/G3_1A.ecc_caller_out.details.nolowq.txt

bedtools intersect -f 1 -wb -a ${GENE_BEDFILE} -b ${ECCDNA_FILE} | awk '{print $5,$6,$7}' | sort | uniq -c | wc -l

MAPFILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile

while read SAMPLE; do
    total=$(bedtools intersect -f 1 -wb -a ${GENE_BEDFILE} -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${SAMPLE}/${SAMPLE}.ecc_caller_out.details.nolowq.txt | awk '{print $5,$6,$7}' | sort | uniq -c | wc -l)
    two_genes=$(bedtools intersect -f 1 -wb -a ${GENE_BEDFILE} -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${SAMPLE}/${SAMPLE}.ecc_caller_out.details.nolowq.txt | awk '{print $5,$6,$7}' | sort | uniq -c | awk '$1 > 1' | wc -l)
    awk -v total=$total -v two=$two_genes '{ print two/total }'
done < ${MAPFILE}