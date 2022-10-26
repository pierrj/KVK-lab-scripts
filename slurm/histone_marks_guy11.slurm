#!/bin/bash
#SBATCH --job-name=histone_marks_guy11
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

# cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/

# ## download and map all chpseq data from Zhang et al 2021
# while read SRA; do
#     /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SRA} -O .
#     /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SRA}.sra
# done < accessions_histone_marks.txt 

# while read SRA; do
#     /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SRA} -O .
#     /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SRA}.sra
# done < chpseq_loading_control_accessions.txt

# # reformat metadata downloaded from Zhang et al
# # sort accessions by treatment
# awk -v FS=',' -v OFS='\t' '{print $1,$50, $55}' metadata_histone_marks.txt > metadata_biorep_treatment.txt
# awk -v OFS='\t' '{print $2, $3}' metadata_biorep_treatment.txt | tail -n+2 | sort | uniq > biorep_treament
# while read biorep treatment; do
# awk -v b=$biorep -v t=$treatment '{if ($2==b && $3==t) {print "raw_data/"$1".sra.fastq"}}' metadata_biorep_treatment.txt > accessions
# cat $(cat accessions | tr "\n" " ") > samples_merged/${treatment}_${biorep}.fastq
# done < biorep_treament

# awk -v FS=',' -v OFS='\t' '{print $1,$50, $55}' chpseq_loading_control_metadata.txt > control_biorep_treatment.txt
# awk -v OFS='\t' '{print $2, $3}' control_biorep_treatment.txt | tail -n+2 | sort | uniq > control_biorep_treament
# while read biorep treatment; do
# awk -v b=$biorep -v t=$treatment '{if ($2==b && $3==t) {print "raw_data/"$1".sra.fastq"}}' control_biorep_treatment.txt > accessions
# cat $(cat accessions | tr "\n" " ") > samples_merged/${treatment}_${biorep}.fastq
# done < control_biorep_treament

# ## map and sort all reads
# while read sample; do
#     bwa mem -t ${SLURM_NTASKS} /global/scratch/users/pierrj/eccDNA/pipeline_tests/eccdna_gene_expression/zhang_et_al_data/guy11_genome_baoetal2017.fasta ${sample}.fastq -o ${sample}.sam
# done < mapfile

# while read sample; do
#     bwa mem -t ${SLURM_NTASKS} /global/scratch/users/pierrj/eccDNA/pipeline_tests/eccdna_gene_expression/zhang_et_al_data/guy11_genome_baoetal2017.fasta ${sample}.fastq -o ${sample}.sam
# done < mapfile_input

# while read biorep treatment; do
#     samtools sort ${treatment}_${biorep}.sam.bam > sorted.${treatment}_${biorep}.sam.bam
#     samtools index sorted.${treatment}_${biorep}.sam.bam
# done < biorep_treament

# while read biorep treatment; do
#     samtools sort ${treatment}_${biorep}.sam.bam > sorted.${treatment}_${biorep}.sam.bam
#     samtools index sorted.${treatment}_${biorep}.sam.bam
# done < control_biorep_treament

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11

## fix guy11 gff to match input for orthofinder

if [ -d GUY1 ]; then
    rm -r GUY11
fi

# mkdir GUY11/
# mkdir GUY11/fungap_out/
# mkdir GUY11/fungap_out/fungap_out/

cp /global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.gff3 GUY11/fungap_out/fungap_out/fungap_out.gff3

MAPFILE=guy11.txt
genome=GUY11

echo GUY11 > $MAPFILE

grep -v '#' ${genome}/fungap_out/fungap_out/fungap_out.gff3 > ${genome}/fungap_out/fungap_out/fungap_out.nocomments.gff3
sed -i '1s/^/##gff-version 3\n/' ${genome}/fungap_out/fungap_out/fungap_out.nocomments.gff3

/global/scratch/users/pierrj/conda_envs/orthofinder/bin/python /global/home/users/pierrj/git/python/process_gffs_for_orthofinder.py $MAPFILE $LINEAGE_INFO

mv all_gffs_fixed/${genome}_fungap_out.fixed.gff3 all_gffs_fixed/${genome}_fungap_out.fixed.gff3.old
grep -v gff-version=3=true all_gffs_fixed/${genome}_fungap_out.fixed.gff3.old > all_gffs_fixed/${genome}_fungap_out.fixed.gff3
rm all_gffs_fixed/${genome}_fungap_out.fixed.gff3.old

cp all_gffs_fixed/${genome}_fungap_out.fixed.gff3 .

genome=$(basename ${genome}_fungap_out.fixed.gff3 | awk '{print substr($1, 0, length($1)-22)}')
awk -v OFS='\t' '{if ($3 == "gene") {print $1, $4, $5, $9}}' ${genome}_fungap_out.fixed.gff3 > ${genome}.genes.bed

awk '{print $4}' ${genome}.genes.bed > ${genome}.genes_column

awk '{print $3-$2}' ${genome}.genes.bed > ${genome}.gene_lengths

cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/biorep_treament .
awk '{print $2}' biorep_treament | sort | uniq > treatments

## run bedtools coverage on the mapped histone mark files, normalize by read count
while read biorep treatment; do
    bam_file=/global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/sorted.${treatment}_${biorep}.sam.bam
    bedtools coverage -counts -a ${genome}.genes.bed -b $bam_file | awk '{print $NF}' > ${treatment}_${biorep}.coverage
    paste ${treatment}_${biorep}.coverage ${genome}.gene_lengths | awk '{print $1/($2)}' > ${treatment}_${biorep}.coverage.RPK
    sum_rpks=$(awk '{sum+=$1;} END{print sum/1000000;}' ${treatment}_${biorep}.coverage.RPK)
    awk -v N=$sum_rpks '{print $1/N}' ${treatment}_${biorep}.coverage.RPK > ${treatment}_${biorep}.coverage.normalized
done < biorep_treament

while read treatment; do
    paste ${treatment}_*.coverage.normalized | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > ${treatment}.coverage.normalized
    paste ${genome}.genes_column ${treatment}.coverage.normalized > guy11_${treatment}_per_gene.txt
done < treatments

## get peak locations and reformat, downloaded from Zhang et al 2021
cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/methyl_peaks_zhang_etal.txt .
cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/acetylation_peaks_zhang_etal.txt .
head -n-1 methyl_peaks_zhang_etal.txt | tail -n+3 | awk -v OFS='\t' '{print $1, $2, $3}' > methyl_peaks_zhang_etal.bed
head -n-1 acetylation_peaks_zhang_etal.txt | tail -n+3 | awk -v OFS='\t' '{print $1, $2, $3}' > acetylation_peaks_zhang_etal.bed

GENOME=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta

# rename files to work with my inputs
grep '>' ${GENOME} | awk -v OFS='\t' '{print substr($1, 2), substr($6,7, length($6)-7)}' > rename_mapfile
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' rename_mapfile methyl_peaks_zhang_etal.bed > methyl_peaks_zhang_etal.renamed.bed
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' rename_mapfile acetylation_peaks_zhang_etal.bed > acetylation_peaks_zhang_etal.renamed.bed

# group genes according to how Zhang et al 2021 dit it
bedtools intersect -wa -f 0.5 -a ${genome}.genes.bed -b acetylation_peaks_zhang_etal.renamed.bed > acetylated_genes
bedtools intersect -v -wa -f 0.5 -a ${genome}.genes.bed -b acetylation_peaks_zhang_etal.renamed.bed > nonacetylated_genes
bedtools intersect -wa -f 0.5 -a ${genome}.genes.bed -b methyl_peaks_zhang_etal.renamed.bed > methylated_genes
bedtools intersect -v -wa -f 0.5 -a ${genome}.genes.bed -b methyl_peaks_zhang_etal.renamed.bed > nonmethylated_genes

cat <(awk -v OFS='\t' '{print $(NF), "TRUE"}' acetylated_genes) <(awk -v OFS='\t' '{print $(NF), "FALSE"}' nonacetylated_genes) > guy11_gene_zhang_acetylation.txt

cat <(awk -v OFS='\t' '{print $(NF), "TRUE"}' methylated_genes) <(awk -v OFS='\t' '{print $(NF), "FALSE"}' nonmethylated_genes) > guy11_gene_zhang_methylation.txt