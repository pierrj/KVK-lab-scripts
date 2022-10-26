#!/bin/bash
#SBATCH --job-name=bismark
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/bismark/moryzae_test

SAMPLE=SRR653493

module load bowtie2
module load samtools

/global/scratch/users/pierrj/bismark/Bismark-0.24.0/bismark_genome_preparation \
    --verbose --bowtie2 --path_to_aligner /global/home/groups/consultsw/sl-7.x86_64/modules/bowtie2/2.3.4.1 guy11_genome_baoetal2017/


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
    --ignore 2 \
    --ignore_r2 2 \
    ${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.bam

CHROM_SIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes

source activate /global/scratch/users/pierrj/conda_envs/deeptools

gunzip ${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.bedGraph.gz

grep -v track ${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.bedGraph > ${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.fixed.bg

bedGraphToBigWig ${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.fixed.bg $CHROM_SIZES ${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.fixed.bw

density_file=${SAMPLE}_bismark_out/${SAMPLE}.sra_1_bismark_bt2_pe.fixed.bw
THREADS=$SLURM_NTASKS
OUTPUT_NAME=methylation
REGIONS_BED=/global/scratch/users/pierrj/PAV_SV/SV/dels_re_rice_blast/all_DEL.bed

computeMatrix scale-regions -p ${THREADS} -S $density_file \
                            -R ${REGIONS_BED} \
                            --beforeRegionStartLength 500 \
                            --regionBodyLength 10 \
                            --afterRegionStartLength 500 \
                            -o ${OUTPUT_NAME}.mat.gz

plotProfile -m ${OUTPUT_NAME}.mat.gz \
            -out ${OUTPUT_NAME}.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "${OUTPUT_NAME}" \
            --outFileNameData ${OUTPUT_NAME}.tab