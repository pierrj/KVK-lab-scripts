#!/bin/bash
#SBATCH --job-name=histone_marks_plotprofile_final
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


## download and map all chpseq data from Zhang et al 2021


# cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/

# while read SRA; do
#     /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SRA} -O .
#     /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SRA}.sra
# done < accessions_histone_marks.txt 

# while read SRA; do
#     /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SRA} -O .
#     /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SRA}.sra
# done < chpseq_loading_control_accessions.txt

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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/histone_marks

## make tracks

# if [ -f "large_eccdnas_all" ]; then
#     rm large_eccdnas_all
# fi

# while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/large_eccdnas >> large_eccdnas_all; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
# sort -k1,1 -k2,2n large_eccdnas_all | uniq > large_eccdnas_all_uniq
# mv large_eccdnas_all_uniq large_eccdnas

# if [ -f "micro_dnas_all" ]; then
#     rm micro_dnas_all
# fi

# while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/micro_dnas >> micro_dnas_all; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
# sort -k1,1 -k2,2n micro_dnas_all | uniq > micro_dnas_all_uniq
# mv micro_dnas_all_uniq micro_dnas

# if [ -f "ltr_eccdnas_all" ]; then
#     rm ltr_eccdnas_all
# fi

# while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/ltr_eccs >> ltr_eccdnas_all; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
# sort -k1,1 -k2,2n ltr_eccdnas_all | uniq > ltr_eccdnas_all_uniq
# mv ltr_eccdnas_all ltr_eccdnas

GENOME=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
GENOME_CHROMSIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes

# cat ltr_eccdnas large_eccdnas micro_dnas > all_eccdnas
# bedtools shuffle -i all_eccdnas -g ${GENOME_CHROMSIZES} > genome_baseline

# cp /global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_high_conf_ltr_tes_for_plotprofile.bed ltr_tes

# cp /global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.bed all_genes


# cp /global/scratch/pierrj/eccDNA/pipeline_tests/histone_marks/methyl_peaks_zhang_etal.txt .
# cp /global/scratch/pierrj/eccDNA/pipeline_tests/histone_marks/acetylation_peaks_zhang_etal.txt .
# head -n-1 methyl_peaks_zhang_etal.txt | tail -n+3 | awk -v OFS='\t' '{print $1, $2, $3}' > methyl_peaks_zhang_etal.bed
# head -n-1 acetylation_peaks_zhang_etal.txt | tail -n+3 | awk -v OFS='\t' '{print $1, $2, $3}' > acetylation_peaks_zhang_etal.bed
# grep '>' ${GENOME} | awk -v OFS='\t' '{print substr($1, 2), substr($6,7, length($6)-7)}' > rename_mapfile

# awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' rename_mapfile methyl_peaks_zhang_etal.bed > methyl_peaks_zhang_etal.renamed.bed
# awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' rename_mapfile acetylation_peaks_zhang_etal.bed > acetylation_peaks_zhang_etal.renamed.bed

# bedtools intersect -wa -f 0.5 -a /global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.bed -b acetylation_peaks_zhang_etal.renamed.bed > acetylated_genes
# bedtools intersect -v -wa -f 0.5 -a /global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.bed -b acetylation_peaks_zhang_etal.renamed.bed > nonacetylated_genes
# bedtools intersect -wa -f 0.5 -a /global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.bed -b methyl_peaks_zhang_etal.renamed.bed > methylated_genes
# bedtools intersect -v -wa -f 0.5 -a /global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.bed -b methyl_peaks_zhang_etal.renamed.bed > nonmethylated_genes

source activate deeptools

data_dir=/global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/

biorep=1

treatment=H3K27ac
# bamCompare -p ${SLURM_NTASKS} -b1 ${data_dir}/sorted.${treatment}_${biorep}.sam.bam -b2 ${data_dir}/sorted.Input_${biorep}.sam.bam -o log2ratio_${treatment}_${biorep}.bw -of bigwig --scaleFactorsMethod readCount
# computeMatrix scale-regions -p ${SLURM_NTASKS} -S log2ratio_${treatment}_${biorep}.bw \
#                             -R genome_baseline all_genes ltr_tes large_eccdnas micro_dnas ltr_eccdnas acetylated_genes nonacetylated_genes \
#                             --beforeRegionStartLength 300 \
#                             --regionBodyLength 1000 \
#                             --afterRegionStartLength 300 \
#                             --skipZeros -o ${treatment}_${biorep}.mat.gz 
plotProfile -m ${treatment}_${biorep}.mat.gz \
            -out ${treatment}_${biorep}.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "${treatment} ${biorep} eccDNAs" \
            --outFileNameData ${treatment}_${biorep}.tab

treatment=H3K27me3
# bamCompare -p ${SLURM_NTASKS} -b1 ${data_dir}/sorted.${treatment}_${biorep}.sam.bam -b2 ${data_dir}/sorted.Input_${biorep}.sam.bam -o log2ratio_${treatment}_${biorep}.bw -of bigwig --scaleFactorsMethod readCount
# computeMatrix scale-regions -p ${SLURM_NTASKS} -S log2ratio_${treatment}_${biorep}.bw \
#                             -R  genome_baseline all_genes ltr_tes large_eccdnas micro_dnas ltr_eccdnas methylated_genes nonmethylated_genes \
#                             --beforeRegionStartLength 300 \
#                             --regionBodyLength 1000 \
#                             --afterRegionStartLength 300 \
#                             --skipZeros -o ${treatment}_${biorep}.mat.gz 
plotProfile -m ${treatment}_${biorep}.mat.gz \
            -out ${treatment}_${biorep}.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "${treatment} ${biorep} eccDNAs" \
            --outFileNameData ${treatment}_${biorep}.tab

treatment=H3K36me3
# bamCompare -p ${SLURM_NTASKS} -b1 ${data_dir}/sorted.${treatment}_${biorep}.sam.bam -b2 ${data_dir}/sorted.Input_${biorep}.sam.bam -o log2ratio_${treatment}_${biorep}.bw -of bigwig --scaleFactorsMethod readCount
# computeMatrix scale-regions -p ${SLURM_NTASKS} -S log2ratio_${treatment}_${biorep}.bw \
#                             -R  genome_baseline all_genes ltr_tes large_eccdnas micro_dnas ltr_eccdnas \
#                             --beforeRegionStartLength 300 \
#                             --regionBodyLength 1000 \
#                             --afterRegionStartLength 300 \
#                             --skipZeros -o ${treatment}_${biorep}.mat.gz 
plotProfile -m ${treatment}_${biorep}.mat.gz \
            -out ${treatment}_${biorep}.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "${treatment} ${biorep} eccDNAs" \
            --outFileNameData ${treatment}_${biorep}.tab

# bedtools makewindows -g ${GENOME_CHROMSIZES} -w 50 > guy11_genome_baoetal2017.50bpwindows

# bedtools nuc -fi ${GENOME} -bed guy11_genome_baoetal2017.50bpwindows > guy11_genome_baoetal2017.50bpwindows.nuc_out

# awk -v OFS='\t' '{ if (NR > 1) {print $1, $2, $3, $5}}' guy11_genome_baoetal2017.50bpwindows.nuc_out > guy11_genome_baoetal2017.50bpwindows.bg

# bedGraphToBigWig guy11_genome_baoetal2017.50bpwindows.bg ${GENOME_CHROMSIZES} guy11_genome_baoetal2017.50bpwindows.bw

# computeMatrix scale-regions -p ${SLURM_NTASKS} -S guy11_genome_baoetal2017.50bpwindows.bw \
#                             -R genome_baseline all_genes ltr_tes large_eccdnas micro_dnas ltr_eccdnas \
#                             --beforeRegionStartLength 300 \
#                             --regionBodyLength 1000 \
#                             --afterRegionStartLength 300 \
#                             --skipZeros -o gc.mat.gz
plotProfile -m gc.mat.gz \
            -out gc.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "GC eccDNAs" \
            --outFileNameData gc.tab