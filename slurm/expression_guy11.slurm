#!/bin/bash
#SBATCH --job-name=expression_guy11
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

THREADS=24

GENOME_FASTA=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
GFF_FILE=GUY11_fungap_out.fixed.gff3
SRA_LIST=zhang_et_al_2019_complete_medium
genome=GUY11

if [ -f "${SRA_LIST}" ]; then
    rm "${SRA_LIST}"
fi

echo SRR12319289 >> $SRA_LIST
echo SRR12319290 >> $SRA_LIST
echo SRR12319292 >> $SRA_LIST
echo SRR12319293 >> $SRA_LIST
echo SRR12319294 >> $SRA_LIST
echo SRR12319295 >> $SRA_LIST
echo SRR12319296 >> $SRA_LIST
echo SRR12319297 >> $SRA_LIST
echo SRR12319298 >> $SRA_LIST
echo SRR12319299 >> $SRA_LIST
echo SRR12319300 >> $SRA_LIST
echo SRR12319301 >> $SRA_LIST

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


GENOME_FASTA=guy11_plus_orsa_w_organelles.fasta
GFF_FILE=GUY11_fungap_out.fixed.w_orsa.gff3
SRA_LIST=zhang_et_al_2019_in_planta
genome=GUY11

if [ -f "${SRA_LIST}" ]; then
    rm "${SRA_LIST}"
fi

echo SRR12319276 >> $SRA_LIST
echo SRR12319277 >> $SRA_LIST
echo SRR12319278 >> $SRA_LIST
echo SRR12319279 >> $SRA_LIST
echo SRR12319281 >> $SRA_LIST
echo SRR12319282 >> $SRA_LIST
echo SRR12319283 >> $SRA_LIST
echo SRR12319284 >> $SRA_LIST
echo SRR12319285 >> $SRA_LIST
echo SRR12319287 >> $SRA_LIST
echo SRR12319288 >> $SRA_LIST
echo SRR12319286 >> $SRA_LIST

basename_gff_file=$(basename ${GFF_FILE})


grep 'MQOP010000' ${GFF_FILE} | grep "exon" | awk -F'[=\t]' -v OFS='\t' '{print $(NF), $5-$4}' | \
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
    awk 'NR>4 {print $0}' ${SRA}.ReadsPerGene.out.tab | grep GUY11 | sort -k1,1 > ${SRA}.ReadsPerGene.out.sorted.tab
    awk '{print $2}' ${SRA}.ReadsPerGene.out.sorted.tab > ${SRA}.ReadsPerGene.out.genecolumn.tab
    paste ${SRA}.ReadsPerGene.out.genecolumn.tab ${basename_gff_file}.exon_lengths | awk '{print $1/($2)}' > ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab
    sum_rpks=$(awk '{sum+=$1;} END{print sum/1000000;}' ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab)
    awk -v N=$sum_rpks '{print $1/N}' ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab > ${SRA}.${SRA_LIST}.TPM.ReadsPerGene.out.genecolumn.tab
done < ${SRA_LIST}

# watch out for wild cards here
paste *.${SRA_LIST}.TPM.ReadsPerGene.out.genecolumn.tab | \
    awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > ${SRA_LIST}.all.TPM.ReadsPerGene.out.genecolumn.tab

paste ${genome}.genes_column ${SRA_LIST}.all.TPM.ReadsPerGene.out.genecolumn.tab > ${genome}_${SRA_LIST}_expression.txt