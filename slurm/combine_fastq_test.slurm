#!/bin/bash
#SBATCH --job-name=combine_fastq_test
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11

THREADS=$SLURM_NTASKS
GENOME_FASTA=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
GFF_FILE=GUY11_fungap_out.fixed.gff3
SRA_LIST=zhang_et_al_2019_complete_medium
genome=GUY11
genome_fasta_basename=$(basename ${GENOME_FASTA})
basename_gff_file=$(basename ${GFF_FILE})

cat SRR12319298.sra.fastq SRR12319299.sra.fastq SRR12319300.sra.fastq SRR12319301.sra.fastq > combined.sra.fastq

SRA=combined

module load STAR

STAR --runThreadN ${THREADS} \
    --genomeDir ${genome_fasta_basename}_starindex \
    --readFilesIn ${SRA}.sra.fastq \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${SRA}. \
    --quantMode GeneCounts
bam_file=${SRA}.Aligned.sortedByCoord.out.bam
awk 'NR>4 {print $0}' ${SRA}.ReadsPerGene.out.tab | sort -k1,1 > ${SRA}.ReadsPerGene.out.sorted.tab
awk '{print $2}' ${SRA}.ReadsPerGene.out.sorted.tab > ${SRA}.ReadsPerGene.out.genecolumn.tab
paste ${SRA}.ReadsPerGene.out.genecolumn.tab ${basename_gff_file}.exon_lengths | awk '{print $1/($2)}' > ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab
sum_rpks=$(awk '{sum+=$1;} END{print sum/1000000;}' ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab)
awk -v N=$sum_rpks '{print $1/N}' ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab > ${SRA}.${SRA_LIST}.TPM.ReadsPerGene.out.genecolumn.tab

FASTA_FILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
MAPFILE=${FASTA_FILE}.contignames

samtools sort -@ ${SLURM_NTASKS} -o sorted.${SRA}.Aligned.sortedByCoord.out.bam ${SRA}.Aligned.sortedByCoord.out.bam
samtools index sorted.${SRA}.Aligned.sortedByCoord.out.bam
samtools view -b sorted.${SRA}.Aligned.sortedByCoord.out.bam $(cat ${MAPFILE} | tr "\n" " ") > filtered.${SRA}.Aligned.sortedByCoord.out.bam

source activate /global/scratch/users/pierrj/conda_envs/deeptools

sv=DEL
FASTA_FILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
GENOME_FILE=$FASTA_FILE
GENOME_CHROMSIZES=${FASTA_FILE}.chromsizes
THREADS=${SLURM_NTASKS}
WINDOWS=10
density_file=filtered.${SRA}.Aligned.sortedByCoord.out.bam
genome_basename=guy11
OUTPUT_NAME=combined
CHROM_SIZES=$GENOME_CHROMSIZES

bedtools makewindows -g ${CHROM_SIZES} -w ${WINDOWS} > ${genome_basename}.${WINDOWS}windows

basename_density_file=$(basename $density_file)
read_count=$(samtools view -c -F 4 -F 2048 $density_file | awk '{print $1/1000000}')
bedtools coverage -counts -sorted -a ${genome_basename}.${WINDOWS}windows \
    -b $density_file -g ${CHROM_SIZES} | awk -v r=$read_count -v OFS='\t' '{print $(NF)/r}' > ${basename_density_file}.${OUTPUT_NAME}.bg

paste *.${OUTPUT_NAME}.bg | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > averaged.${OUTPUT_NAME}.bg
paste ${genome_basename}.${WINDOWS}windows averaged.${OUTPUT_NAME}.bg > ${density_file_basename}.bg
bedGraphToBigWig ${density_file_basename}.bg ${CHROM_SIZES} ${density_file_basename}.bw

sv=DEL
REGIONS_BED=/global/scratch/users/pierrj/PAV_SV/SV/dels_re_rice_blast/shuffled_${sv}.bed

computeMatrix scale-regions -p ${THREADS} -S ${density_file_basename}.bw \
                            -R ${REGIONS_BED} \
                            --beforeRegionStartLength 1000 \
                            --regionBodyLength 500 \
                            --afterRegionStartLength 1000 \
                            -o ${OUTPUT_NAME}.shuffled.mat.gz

sv=DEL
REGIONS_BED=/global/scratch/users/pierrj/PAV_SV/SV/dels_re_rice_blast/all_genes.bed

computeMatrix scale-regions -p ${THREADS} -S ${density_file_basename}.bw \
                            -R ${REGIONS_BED} \
                            --beforeRegionStartLength 1000 \
                            --regionBodyLength 500 \
                            --afterRegionStartLength 1000 \
                            -o ${OUTPUT_NAME}.genes.mat.gz