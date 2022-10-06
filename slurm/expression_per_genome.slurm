#!/bin/bash
#SBATCH --job-name=expression_per_genome
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


## launch with
# genome=something
# sbatch --export=ALL,genome=$genome /global/home/users/pierrj/git/slurm/expression_per_genome.slurm


PROJ_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/${genome}_expression
THREADS=$SLURM_NTASKS

if [ -d "$PROJ_DIR" ]; then
    rm -r "$PROJ_DIR"
fi

mkdir $PROJ_DIR

cd $PROJ_DIR

GENOME_FASTA=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/assemblies/${genome}.fasta
GFF_FILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/all_gffs_fixed/${genome}_fungap_out.fixed.gff3
SRA_LIST=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11/zhang_et_al_2019_complete_medium
SRA_TAG=$(basename $SRA_LIST)

grep 'exon' ${GFF_FILE} | awk -F'[=\t]' -v OFS='\t' '{print $(NF), $5-$4}' | \
    awk '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }' | \
    sort -k1,1 | awk '{print $1}' > ${genome}.genes_column

basename_gff_file=$(basename ${GFF_FILE})

grep 'exon' ${GFF_FILE} | awk -F'[=\t]' -v OFS='\t' '{print $(NF), $5-$4}' | \
    awk '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }' | \
    sort -k1,1 | awk '{print $2/1000}' > ${basename_gff_file}.exon_lengths

genome_fasta_basename=$(basename ${GENOME_FASTA})

if [ -d "${genome_fasta_basename}_starindex" ]; then
    rm -r ${genome_fasta_basename}_starindex
fi

mkdir ${genome_fasta_basename}_starindex

module load STAR

## index genome for STAR
STAR --runThreadN ${THREADS} --runMode genomeGenerate --genomeDir ${genome_fasta_basename}_starindex \
    --genomeFastaFiles ${GENOME_FASTA} \
    --sjdbGTFfile ${GFF_FILE} \
    --sjdbOverhang 100 \
    --genomeSAindexNbases 11 \
    --sjdbGTFtagExonParentTranscript ID \
    --sjdbGTFtagExonParentGene Parent

while read SRA; do
    STAR --runThreadN ${THREADS} \
        --genomeDir ${genome_fasta_basename}_starindex \
        --readFilesIn /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11/${SRA}.sra.fastq \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${SRA}. \
        --quantMode GeneCounts
    bam_file=${SRA}.Aligned.sortedByCoord.out.bam
    awk 'NR>4 {print $0}' ${SRA}.ReadsPerGene.out.tab | sort -k1,1 > ${SRA}.ReadsPerGene.out.sorted.tab
    awk '{print $2}' ${SRA}.ReadsPerGene.out.sorted.tab > ${SRA}.ReadsPerGene.out.genecolumn.tab
    num_reads=$(samtools view -c -F 4 -F 2048 $bam_file | awk '{print $1/1000000}')
    paste ${SRA}.ReadsPerGene.out.genecolumn.tab ${basename_gff_file}.exon_lengths | awk -v N=$num_reads '{print $1/($2*N)}' > ${SRA}.${SRA_TAG}.RPKM.ReadsPerGene.out.genecolumn.tab
done < ${SRA_LIST}

paste *.${SRA_TAG}.RPKM.ReadsPerGene.out.genecolumn.tab | \
    awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > all.${SRA_TAG}.RPKM.ReadsPerGene.out.genecolumn.tab

paste ${genome}.genes_column all.${SRA_TAG}.RPKM.ReadsPerGene.out.genecolumn.tab > ${genome}_${SRA_TAG}_expression.txt



## FOR IN PLANTA EXPRESSION ##

RICE_GENOME_FASTA=/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0.fasta
cat $RICE_GENOME_FASTA $GENOME_FASTA > ${genome}_plus_rice.fasta

RICE_GFF=/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0.gff
cat $RICE_GFF $GFF_FILE > ${genome}_plus_rice.gff

GFF_FILE=${genome}_plus_rice.gff
GENOME_FILE=${genome}_plus_rice.fasta
SRA_LIST=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11/zhang_et_al_2019_in_planta
SRA_TAG=$(basename $SRA_LIST)

basename_gff_file=$(basename ${GFF_FILE})


grep 'Scaffold' ${GFF_FILE} | grep "exon" | awk -F'[=\t]' -v OFS='\t' '{print $(NF), $5-$4}' | \
    awk '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }' | \
    sort -k1,1 | awk '{print $2/1000}' > ${basename_gff_file}.exon_lengths

genome_fasta_basename=$(basename ${GENOME_FASTA})

if [ -d "${genome_fasta_basename}_starindex" ]; then
    rm -r ${genome_fasta_basename}_starindex
fi

mkdir ${genome_fasta_basename}_starindex

module load STAR

## index genome for STAR
STAR --runThreadN ${THREADS} --runMode genomeGenerate --genomeDir ${genome_fasta_basename}_starindex \
    --genomeFastaFiles ${GENOME_FASTA} \
    --sjdbGTFfile ${GFF_FILE} \
    --sjdbOverhang 100 \
    --genomeSAindexNbases 11 \
    --sjdbGTFtagExonParentTranscript ID \
    --sjdbGTFtagExonParentGene Parent

while read SRA; do
    STAR --runThreadN ${THREADS} \
        --genomeDir ${genome_fasta_basename}_starindex \
        --readFilesIn /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11/${SRA}.sra.fastq \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${SRA}. \
        --quantMode GeneCounts
    bam_file=${SRA}.Aligned.sortedByCoord.out.bam
    awk 'NR>4 {print $0}' ${SRA}.ReadsPerGene.out.tab | grep ${genome} | sort -k1,1 > ${SRA}.ReadsPerGene.out.sorted.tab
    awk '{print $2}' ${SRA}.ReadsPerGene.out.sorted.tab > ${SRA}.ReadsPerGene.out.genecolumn.tab
    num_reads=$(samtools view -c -F 4 -F 2048 $bam_file | awk '{print $1/1000000}')
    paste ${SRA}.ReadsPerGene.out.genecolumn.tab ${basename_gff_file}.exon_lengths | awk -v N=$num_reads '{print $1/($2*N)}' > ${SRA}.${SRA_TAG}.RPKM.ReadsPerGene.out.genecolumn.tab
done < ${SRA_LIST}

paste *.${SRA_TAG}.RPKM.ReadsPerGene.out.genecolumn.tab | \
    awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > all.${SRA_TAG}.RPKM.ReadsPerGene.out.genecolumn.tab

paste ${genome}.genes_column all.${SRA_TAG}.RPKM.ReadsPerGene.out.genecolumn.tab > ${genome}_${SRA_TAG}_expression.txt