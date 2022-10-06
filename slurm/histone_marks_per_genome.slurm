#!/bin/bash
#SBATCH --job-name=histone_marks_per_genome
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

## launch with
# genome=something
# sbatch --export=ALL,genome=$genome /global/home/users/pierrj/git/slurm/expression_per_genome.slurm

PROJ_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/${genome}_histone_marks
THREADS=$SLURM_NTASKS

if [ -d "$PROJ_DIR" ]; then
    rm -r "$PROJ_DIR"
fi

mkdir $PROJ_DIR

cd $PROJ_DIR

GENOME_FASTA=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/assemblies/${genome}.fasta
GFF_FILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/all_gffs_fixed/${genome}_fungap_out.fixed.gff3
MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/mapfile
TREATMENT_MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/mapfile

awk -v OFS='\t' '{if ($3 == "gene") {print $1, $4, $5, $9}}' $GFF_FILE > ${genome}.genes.bed
awk '{print $4}' ${genome}.genes.bed > ${genome}.genes_column


## map and sort all reads
while read sample; do
    bwa mem -t ${THREADS} $GENOME_FASTA /global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/${sample}.fastq | \
        samtools sort -@ ${THREADS} -o ${sample}.bam -
    samtools index ${sample}.bam
    read_count=$(samtools view -c -F 4 -F 2048 ${sample}.bam | awk '{print $1}')
    ## run bedtools coverage on the mapped histone mark files, normalize by read count
    bedtools coverage -counts -a ${genome}.genes.bed -b $bam_file | awk -v r=$read_count '{print $(NF)/(r*1000000)}' > ${sample}.coverage
done < $MAPFILE

while read treatment; do
    paste ${treatment}_*.coverage | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > ${treatment}.coverage
    paste ${genome}.genes_column ${treatment}.coverage > ${genome}_${treatment}_per_gene.txt
done < $TREATMENT_MAPFILE