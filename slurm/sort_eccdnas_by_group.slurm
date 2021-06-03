#!/bin/bash
#SBATCH --job-name=sort_eccdnas_by_group
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina_if_v_g3

gff_file=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.gff3
copia_file=/global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_copia_locs.bed
gypsy_file=/global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_gypsy_locs.bed

while read sample; do
cd ${sample}
        awk '{if ($3 == "gene") print $0}' ${gff_file} > tmp.genes.gff
        awk ' $5 != "lowq" ' ${sample}.ecc_caller_out.details.txt > ${sample}.ecc_caller_out.details.nolowq.txt
        bedtools intersect -u -f 0.5 -a ${sample}.ecc_caller_out.details.nolowq.txt -b tmp.genes.gff | awk -v OFS='\t' '{print $1, $2, $3}' > genic_eccs
        sed 's/[[:space:]]*$//' ${sample}.ecc_caller_out.splitreads.bed | awk 'NR==FNR{a[$0]++;next}a[$0]' genic_eccs - > genic_eccs_splitreads
        bedtools intersect -v -f 0.5 -a ${sample}.ecc_caller_out.details.nolowq.txt -b tmp.genes.gff | awk -v OFS='\t' '{print $1, $2, $3}' > noncoding_eccs
        sed 's/[[:space:]]*$//' ${sample}.ecc_caller_out.splitreads.bed | awk 'NR==FNR{a[$0]++;next}a[$0]' noncoding_eccs - > noncoding_eccs_splitreads
        bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b ${copia_file} | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) >= 0.9) {print $1, $2, $3}}' | sort -k1,1 -k2,2n > copia_eccs
        sed 's/[[:space:]]*$//' ${sample}.ecc_caller_out.splitreads.bed | awk 'NR==FNR{a[$0]++;next}a[$0]' copia_eccs - > copia_eccs_splitreads
        bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b ${gypsy_file} | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) >= 0.9) {print $1, $2, $3}}' | sort -k1,1 -k2,2n > gypsy_eccs
        sed 's/[[:space:]]*$//' ${sample}.ecc_caller_out.splitreads.bed | awk 'NR==FNR{a[$0]++;next}a[$0]' gypsy_eccs - > gypsy_eccs_splitreads
        cat copia_eccs gypsy_eccs > ltr_eccs
        sed 's/[[:space:]]*$//' ${sample}.ecc_caller_out.splitreads.bed | awk 'NR==FNR{a[$0]++;next}a[$0]' ltr_eccs - > ltr_eccs_splitreads
        bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b ${copia_file} | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $1, $2, $3}}' | \
            bedtools intersect -wao -a - -b ${gypsy_file} | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $1, $2, $3}}'| sort -k1,1 -k2,2n | awk '$3-$2 <= 400' > micro_dnas
        sed 's/[[:space:]]*$//' ${sample}.ecc_caller_out.splitreads.bed | awk 'NR==FNR{a[$0]++;next}a[$0]' micro_dnas - > micro_dnas_splitreads
        bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b ${copia_file} | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $1, $2, $3}}' | \
            bedtools intersect -wao -a - -b ${gypsy_file} | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $1, $2, $3}}'| sort -k1,1 -k2,2n | awk '$3-$2 > 400' > large_eccdnas
        sed 's/[[:space:]]*$//' ${sample}.ecc_caller_out.splitreads.bed | awk 'NR==FNR{a[$0]++;next}a[$0]' large_eccdnas - > large_eccdnas_splitreads
cd ..
done < mapfile