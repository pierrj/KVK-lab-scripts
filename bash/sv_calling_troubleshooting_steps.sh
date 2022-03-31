

for sv in DEL DUP INV TRA; do
    echo ${sv}
    echo ${sv}
    echo ${sv}
    echo ${sv}
    echo ${sv}
    echo ${sv}
    echo ${sv}
    echo ${sv}
    echo ${sv}
    echo ${sv}
    echo ${sv}
    echo ${sv}
    echo ${sv}    
    awk '$2 > $3' all_${sv}.bed | head -10
done

sv=INV

while read SAMPLE; do
    echo ${SAMPLE}
    awk -F "\t|;" -v OFS='\t' -v sv=${sv} '{ if ( $11 == "SVTYPE="sv) {print $1, $2, substr($14, 5, length($14)-4)}}' \
        /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/${SAMPLE}/all/${SAMPLE}.all.vcf | awk '$2 > $3'
done < /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/gladieux_accessions

awk -F "\t|;" -v OFS='\t' -v sv=${sv} '$11 == "SVTYPE="sv && $1 == "MQOP01000004.1" && $2 == "4338375"' /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/${SAMPLE}/all/${SAMPLE}.all.vcf

MQOP01000004.1  4375532 4374719

## SURVIVOR
sv=DEL

if [ -f all_svs_list ]; then
    rm all_svs_list
fi

while read SAMPLE; do
    echo /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/${SAMPLE}/all/${SAMPLE}.all.vcf >> all_svs_list
done < gladieux_accessions

while read SAMPLE; do
    awk -F "\t|;" -v OFS='\t' -v sv=${sv} '$11 == "SVTYPE="sv' dels_${SAMPLE}.vcf

    echo /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/${SAMPLE}/all/${SAMPLE}.all.vcf >> all_svs_list
done < gladieux_accessions

/global/scratch/users/pierrj/SURVIVOR/Debug/SURVIVOR merge all_svs_list 1000 1 1 1 0 0 all_svs.vcf

awk -F "\t|;" -v OFS='\t' -v sv=${sv} '{ if ( $11 == "SVTYPE="sv) {print $1, $2, substr($14, 5, length($14)-4)}}' all_svs.vcf > all_${sv}.bed

mv all_${sv}.bed all_${sv}.bed.old

grep -v 'Supercontig_7.9' all_${sv}.bed.old > all_${sv}.bed

mv all_${sv}.bed all_${sv}.survivor.bed

# cat
if [ -f all_${sv}.bed ]; then
    rm all_${sv}.bed
fi

while read SAMPLE; do
    awk -F "\t|;" -v OFS='\t' -v sv=${sv} '{ if ( $11 == "SVTYPE="sv) {print $1, $2, substr($14, 5, length($14)-4)}}' \
        /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/${SAMPLE}/all/${SAMPLE}.all.vcf >> all_${sv}.bed
done < gladieux_accessions

mv all_${sv}.bed all_${sv}.bed.old

sort -k1,1 -k2,2n all_${sv}.bed.old | grep -v 'Supercontig_7.9' | uniq > all_${sv}.bed

# sort -k1,1 -k2,2n all_${sv}.bed.old | grep -v 'Supercontig_7.9'  > all_${sv}.bed


mv all_${sv}.bed all_${sv}.cat.bed

###

while read SAMPLE; do
    echo ${SAMPLE}
    awk -F "\t|;" -v OFS='\t' -v sv=${sv} '$11 == "SVTYPE="sv && $1 == "MQOP01000001.1" && $2 == "96186"' /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/${SAMPLE}/all/${SAMPLE}.all.vcf
done < gladieux_accessions


while read SAMPLE; do
    echo ${SAMPLE}
    head -100 /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/${SAMPLE}/all/${SAMPLE}.all.vcf
done < gladieux_accessions

problematic sample: SRR6669187

SAMPLE=SRR6669187
cd /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/${SAMPLE}/

## with strandedness enforced
/global/scratch/users/pierrj/SURVIVOR/Debug/SURVIVOR merge ${SAMPLE}/all/${SAMPLE}.survivor.vcflist 1000 3 1 1 0 50 ${SAMPLE}/all/test.all.vcf

## without strandedness enforced
/global/scratch/users/pierrj/SURVIVOR/Debug/SURVIVOR merge ${SAMPLE}/all/${SAMPLE}.survivor.vcflist 1000 3 1 0 0 50 ${SAMPLE}/all/${SAMPLE}.all.vcf

## vcf2bed tests
/global/scratch/users/pierrj/SURVIVOR/Debug/SURVIVOR vcftobed all_svs.vcf 1 10000000000000000 all_svs.bed

awk '$11 == "DEL"' all_svs.bed > all_DEL.survivor.bed

sv=DEL
awk -F "\t|;" -v OFS='\t' -v sv=${sv} '{ if ( $11 == "SVTYPE="sv) {print $1, $2, substr($14, 5, length($14)-4)}}' all_svs.vcf > all_DEL.mine.bed