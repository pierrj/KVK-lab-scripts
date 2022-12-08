#!/bin/bash
#SBATCH --job-name=expression_b71
#SBATCH --partition=savio3_htc
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all/random_forest/b71

THREADS=$SLURM_NTASKS

GENOME_FASTA=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_to_annotate/GCA_004785725.2_ASM478572v2_genomic.fna
GFF_FILE=/global/scratch/users/pierrj/fungap_runs/wheat_blast/all_gffs_fixed/GCA_004785725.2_ASM478572v2_fungap_out.fixed.gff3
SRA_LIST=b71_assembly_complete_medium
genome=B71

if [ -f "${SRA_LIST}" ]; then
    rm "${SRA_LIST}"
fi

echo SRR9127598 >> $SRA_LIST
echo SRR9127597 >> $SRA_LIST
echo SRR9127600 >> $SRA_LIST

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
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SRA} -O .
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SRA}.sra
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
done < ${SRA_LIST}

# watch out for wild cards here
paste *.${SRA_LIST}.TPM.ReadsPerGene.out.genecolumn.tab | \
    awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > ${SRA_LIST}.all.TPM.ReadsPerGene.out.genecolumn.tab

paste ${genome}.genes_column ${SRA_LIST}.all.TPM.ReadsPerGene.out.genecolumn.tab > ${genome}_${SRA_LIST}_expression.txt

cd /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all/random_forest/b71

THREADS=$SLURM_NTASKS

STAR_INDEX=/global/scratch/projects/vector_kvklab/references/STAR_indices/Triticum_aestivum.IWGSC_plus_B71_starindex/
GFF_FILE=/global/scratch/projects/vector_kvklab/references/gff/chinese_spring_w_b71.gff3
SRA_LIST=b71_assembly_complete_medium_in_planta
genome=B71

if [ -f "${SRA_LIST}" ]; then
    rm "${SRA_LIST}"
fi

echo SRR9127599 >> $SRA_LIST
echo SRR9127602 >> $SRA_LIST
echo SRR9127601 >> $SRA_LIST

basename_gff_file=$(basename ${GFF_FILE})

grep 'CP06033' ${GFF_FILE} | grep "exon" | awk -F'[=\t]' -v OFS='\t' '{print $(NF), $5-$4}' | \
    awk '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }' | \
    sort -k1,1 | awk '{print $2/1000}' > ${basename_gff_file}.exon_lengths

module load STAR

while read SRA; do
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SRA} -O .
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SRA}.sra
    STAR --runThreadN ${THREADS} \
        --genomeDir $STAR_INDEX \
        --readFilesIn ${SRA}.sra.fastq \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${SRA}. \
        --quantMode GeneCounts
    bam_file=${SRA}.Aligned.sortedByCoord.out.bam
    awk 'NR>4 {print $0}' ${SRA}.ReadsPerGene.out.tab | grep GCA004785725.2 | sort -k1,1 > ${SRA}.ReadsPerGene.out.sorted.tab
    awk '{print $2}' ${SRA}.ReadsPerGene.out.sorted.tab > ${SRA}.ReadsPerGene.out.genecolumn.tab
    paste ${SRA}.ReadsPerGene.out.genecolumn.tab ${basename_gff_file}.exon_lengths | awk '{print $1/($2)}' > ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab
    sum_rpks=$(awk '{sum+=$1;} END{print sum/1000000;}' ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab)
    awk -v N=$sum_rpks '{print $1/N}' ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab > ${SRA}.${SRA_LIST}.TPM.ReadsPerGene.out.genecolumn.tab
done < ${SRA_LIST}

# watch out for wild cards here
paste *.${SRA_LIST}.TPM.ReadsPerGene.out.genecolumn.tab | \
    awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > ${SRA_LIST}.all.TPM.ReadsPerGene.out.genecolumn.tab

paste ${genome}.genes_column ${SRA_LIST}.all.TPM.ReadsPerGene.out.genecolumn.tab > ${genome}_${SRA_LIST}_expression.txt