#!/bin/bash
#SBATCH --job-name=methylation_guy11
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


genome=guy11


source activate /global/scratch/users/pierrj/conda_envs/deeptools

module load bowtie2
module load samtools

PROJ_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/${genome}_methylation
THREADS=$SLURM_NTASKS

if [ -d "$PROJ_DIR" ]; then
    rm -r "$PROJ_DIR"
fi

mkdir $PROJ_DIR

cd $PROJ_DIR

GENOME_FASTA=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
GENOME_CHROMSIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes
GFF_FILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11/GUY11_fungap_out.fixed.gff3


mkdir ${genome}_dir

cp $GENOME_FASTA ${genome}_dir/

/global/scratch/users/pierrj/bismark/Bismark-0.24.0/bismark_genome_preparation \
    --verbose --bowtie2 --path_to_aligner /global/home/groups/consultsw/sl-7.x86_64/modules/bowtie2/2.3.4.1 ${genome}_dir/

awk -v OFS='\t' '{if ($3 == "gene") {print $1, $4, $5, $9}}' $GFF_FILE > ${genome}.genes.bed

for sample in SRR653493 SRR653494 SRR653495; do
    /global/scratch/users/pierrj/bismark/Bismark-0.24.0/bismark \
        --parallel 6 \
        --temp_dir tmp \
        --output_dir ${SAMPLE}_bismark_out \
        -p 4 \
        -1 ${SAMPLE}.sra_1.fastq -2 ${SAMPLE}.sra_2.fastq \
        guy11_genome_baoetal2017/

    /global/scratch/users/pierrj/bismark/Bismark-0.24.0/deduplicate_bismark -p \
            --output_dir  ${SAMPLE}_bismark_out \
            ${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.bam

    /global/scratch/users/pierrj/bismark/Bismark-0.24.0//bismark_methylation_extractor -p \
        --output ${SAMPLE}_bismark_out \
        --comprehensive \
        --parallel 10 \
        --bedGraph \
        --CX_context \
        --ignore 2 \
        --ignore_r2 2 \
        ${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.bam

    gunzip ${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.bedGraph.gz

    grep -v track ${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.bedGraph > ${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.fixed.bg

    bedtools intersect -wa -wb -a ${genome}.genes.bed -b ${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.fixed.bg | \
        awk -v OFS='\t' '{a[$4] += $(NF); b[$4] += 1} END{for (i in a) print i, a[i]/b[i]}' > ${genome}.${sample}.percentmethylated

    bedGraphToBigWig ${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.fixed.bg $CHROM_SIZES ${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.fixed.bw

    density_file=${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.fixed.bw
    THREADS=$SLURM_NTASKS
    OUTPUT_NAME=${SAMPLE}_methylation
    REGIONS_BED=/global/scratch/users/pierrj/PAV_SV/SV/dels_re_rice_blast/all_DEL.bed

    computeMatrix scale-regions -p ${THREADS} -S $density_file \
                                -R ${REGIONS_BED} \
                                --beforeRegionStartLength 1000 \
                                --regionBodyLength 500 \
                                --afterRegionStartLength 1000 \
                                -o ${OUTPUT_NAME}.mat.gz

    plotProfile -m ${OUTPUT_NAME}.mat.gz \
                -out ${OUTPUT_NAME}.pdf \
                --numPlotsPerRow 1 \
                --plotTitle "${OUTPUT_NAME}" \
                --outFileNameData ${OUTPUT_NAME}.tab

done