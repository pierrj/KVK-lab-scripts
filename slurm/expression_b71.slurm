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
#MIT License
#
#Copyright (c) 2023 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

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

## get names of all genes
grep 'exon' ${GFF_FILE} | awk -F'[=\t]' -v OFS='\t' '{print $(NF), $5-$4}' | \
    awk '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }' | \
    sort -k1,1 | awk '{print $1}' > ${genome}.genes_column

basename_gff_file=$(basename ${GFF_FILE})

## get lenght of exons
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
    ## download reads
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SRA} -O .
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SRA}.sra
    # map RNAseq reads
    STAR --runThreadN ${THREADS} \
        --genomeDir ${genome_fasta_basename}_starindex \
        --readFilesIn ${SRA}.sra.fastq \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${SRA}. \
        --quantMode GeneCounts
    bam_file=${SRA}.Aligned.sortedByCoord.out.bam
    ## parse reads per gene
    awk 'NR>4 {print $0}' ${SRA}.ReadsPerGene.out.tab | sort -k1,1 > ${SRA}.ReadsPerGene.out.sorted.tab
    awk '{print $2}' ${SRA}.ReadsPerGene.out.sorted.tab > ${SRA}.ReadsPerGene.out.genecolumn.tab
    # calculate RPK values for all genes
    paste ${SRA}.ReadsPerGene.out.genecolumn.tab ${basename_gff_file}.exon_lengths | awk '{print $1/($2)}' > ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab
    # sum up RPKs total
    sum_rpks=$(awk '{sum+=$1;} END{print sum/1000000;}' ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab)
    # calculate per million RPKs
    awk -v N=$sum_rpks '{print $1/N}' ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab > ${SRA}.${SRA_LIST}.TPM.ReadsPerGene.out.genecolumn.tab
done < ${SRA_LIST}

# average all samples together
# watch out for wild cards here
paste *.${SRA_LIST}.TPM.ReadsPerGene.out.genecolumn.tab | \
    awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > ${SRA_LIST}.all.TPM.ReadsPerGene.out.genecolumn.tab

paste ${genome}.genes_column ${SRA_LIST}.all.TPM.ReadsPerGene.out.genecolumn.tab > ${genome}_${SRA_LIST}_expression.txt

cd /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all/random_forest/b71

THREADS=$SLURM_NTASKS

STAR_INDEX=/global/scratch/projects/vector_kvklab/references/STAR_indices/Triticum_aestivum.IWGSC_plus_B71_starindex/
## previously generated star index file
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

## get lenght of exons, dealing with the fact that wheat genes are in there too
grep 'CP06033' ${GFF_FILE} | grep "exon" | awk -F'[=\t]' -v OFS='\t' '{print $(NF), $5-$4}' | \
    awk '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }' | \
    sort -k1,1 | awk '{print $2/1000}' > ${basename_gff_file}.exon_lengths

module load STAR

while read SRA; do
    # download files
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SRA} -O .
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SRA}.sra
    # map reads
    STAR --runThreadN ${THREADS} \
        --genomeDir $STAR_INDEX \
        --readFilesIn ${SRA}.sra.fastq \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${SRA}. \
        --quantMode GeneCounts
    bam_file=${SRA}.Aligned.sortedByCoord.out.bam
    # process gene count file
    awk 'NR>4 {print $0}' ${SRA}.ReadsPerGene.out.tab | grep GCA004785725.2 | sort -k1,1 > ${SRA}.ReadsPerGene.out.sorted.tab
    awk '{print $2}' ${SRA}.ReadsPerGene.out.sorted.tab > ${SRA}.ReadsPerGene.out.genecolumn.tab
    # calculate RPK per gene
    paste ${SRA}.ReadsPerGene.out.genecolumn.tab ${basename_gff_file}.exon_lengths | awk '{print $1/($2)}' > ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab
    # sum up total rpks
    sum_rpks=$(awk '{sum+=$1;} END{print sum/1000000;}' ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab)
    # calculate per million RPKs
    awk -v N=$sum_rpks '{print $1/N}' ${SRA}.${SRA_LIST}.RPK.ReadsPerGene.out.genecolumn.tab > ${SRA}.${SRA_LIST}.TPM.ReadsPerGene.out.genecolumn.tab
done < ${SRA_LIST}

# average all samples together
# watch out for wild cards here
paste *.${SRA_LIST}.TPM.ReadsPerGene.out.genecolumn.tab | \
    awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > ${SRA_LIST}.all.TPM.ReadsPerGene.out.genecolumn.tab

paste ${genome}.genes_column ${SRA_LIST}.all.TPM.ReadsPerGene.out.genecolumn.tab > ${genome}_${SRA_LIST}_expression.txt