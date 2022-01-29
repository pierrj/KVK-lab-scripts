#!/bin/bash
#SBATCH --job-name=plot_profiler_test
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/plot_profiler

repeatmasker/4.1.1

# awk -F "\t|;" -v OFS='\t' '{ if ( $5 == "<DEL>") {print $1, $2, substr($14, 5, length($14)-4)}}' /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/ERR4832458/all/ERR4832458.all.vcf > ERR4832458.dels.bed
# RM_LIB=/global/scratch/users/pierrj/moryzae_pav/te_annotations_anne/anne_repbase_te_lib
# RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . guy11_genome_baoetal2017_with_70-15_mito.fasta

REGIONS_BED=ERR4832458.dels.bed
GENOME_FILE=guy11_genome_baoetal2017_with_70-15_mito.fasta
THREADS=${SLURM_NTASKS}
WINDOWS=10

source activate /global/scratch/users/pierrj/conda_envs/deeptools

DENSITY_BED=guy11_genome_baoetal2017_with_70-15_mito.fasta.out.gff
OUTPUT_NAME=tes_skip_zeroes
/global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o $OUTPUT_NAME


for OUTPUT_NAME in H3K27ac H3K36me3 H3K27me3; do
    TREATMENT_BAM=/global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/sorted.${OUTPUT_NAME}_1.sam.bam
    CONTROL_BAM=/global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/sorted.Input_1.sam.bam
    /global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $TREATMENT_BAM -d $CONTROL_BAM -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_skip_zeroes
done

DENSITY_BED=gc
OUTPUT_NAME=gc_skip_zeroes
/global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o $OUTPUT_NAME