#!/bin/bash
#SBATCH --job-name=plot_profiler_all_svs
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


cd /global/scratch/users/pierrj/PAV_SV/SV/plot_profile

if [ -f all_dels.bed ]; then
    rm all_dels.bed
fi

while read SAMPLE; do
    awk -F "\t|;" -v OFS='\t' '{ if ( $5 == "<DEL>") {print $1, $2, substr($14, 5, length($14)-4)}}' \
        /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/${SAMPLE}/all/${SAMPLE}.all.vcf >> all_dels.bed
done < /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/gladieux_accessions


mv all_dels.bed all_dels.bed.old


## maybe this should be done with survivor instead

sort -k1,1 -k2,2n all_dels.bed.old | grep -v 'Supercontig_7.9' | uniq > all_dels.bed


cp /global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito.fasta .

module load repeatmasker/4.1.1

RM_LIB=/global/scratch/users/pierrj/moryzae_pav/te_annotations_anne/anne_repbase_te_lib
RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . guy11_genome_baoetal2017_with_70-15_mito.fasta

OUTPUT_NAME=plot_profile_all_dels

REGIONS_BED=all_dels.bed
GENOME_FILE=guy11_genome_baoetal2017_with_70-15_mito.fasta
THREADS=${SLURM_NTASKS}
WINDOWS=10

source activate /global/scratch/users/pierrj/conda_envs/deeptools

DENSITY_BED=guy11_genome_baoetal2017_with_70-15_mito.fasta.out.gff
SUB_OUTPUT_NAME=tes
/global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME}

for treatment in H3K27ac H3K36me3 H3K27me3; do
    TREATMENT_BAM=/global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/sorted.${treatment}_1.sam.bam
    CONTROL_BAM=/global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/sorted.Input_1.sam.bam
    /global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $TREATMENT_BAM -d $CONTROL_BAM -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${treatment}
done

DENSITY_BED=gc
SUB_OUTPUT_NAME=gc
/global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME}