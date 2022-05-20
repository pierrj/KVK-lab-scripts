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


# if [ ! -f "guy11_genome_baoetal2017_with_70-15_mito.fasta.out.gff" ]; then
#     module load repeatmasker/4.1.1
#     RM_LIB=/global/scratch/users/pierrj/PAV_SV/PAV/te_annotations_anne/anne_repbase_te_lib_fixed
#     RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . guy11_genome_baoetal2017_with_70-15_mito.fasta
# fi

# if [ -f all_svs_list ]; then
#     rm all_svs_list
# fi
# while read SAMPLE; do
#     echo /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/${SAMPLE}/all/${SAMPLE}.all.vcf >> all_svs_list
# done < /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/gladieux_accessions

# ## in order
# # distance between breakpoints
# # number of supporting callers
# # take type into account
# # take strand into account
# # disabled?
# # minimum size
# # output

# /global/scratch/users/pierrj/SURVIVOR/Debug/SURVIVOR merge all_svs_list 1000 1 1 1 0 0 all_svs.vcf

# /global/scratch/users/pierrj/SURVIVOR/Debug/SURVIVOR vcftobed all_svs.vcf 1 10000000000000000 all_svs.bed

# mv all_svs.bed all_svs.bed.old

# grep -v 'Supercontig_7.9' all_svs.bed.old > all_svs.bed

# for sv in DEL DUP INV TRA; do

#     if [ -f all_${sv}.bed ]; then
#         rm all_${sv}.bed
#     fi

#     awk -v sv=$sv -v OFS='\t' '{ if ( $11 == sv ) {print $12, $13, $15}}' all_svs.bed > all_${sv}.bed

#     # third field for TRAs are all wrong so do breakpoints instead.
#     if [[ "${sv}" == "TRA" ]]
#     then
#         mv all_TRA.bed all_TRA.bed.old
#         awk -v OFS='\t' '{print $1, $2, $2+1}' all_TRA.bed.old > all_TRA.bed
#     fi

#     OUTPUT_NAME=plot_profile_all_${sv}_no_survivor

#     REGIONS_BED=all_${sv}.bed
#     GENOME_FILE=guy11_genome_baoetal2017_with_70-15_mito.fasta
#     THREADS=${SLURM_NTASKS}
#     WINDOWS=10

#     source activate /global/scratch/users/pierrj/conda_envs/deeptools

#     DENSITY_BED=guy11_genome_baoetal2017_with_70-15_mito.fasta.out.gff
#     SUB_OUTPUT_NAME=tes
#     /global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv}

#     DENSITY_BED=guy11_genome_baoetal2017_with_70-15_mito.fasta.simplerepeats.bed
#     SUB_OUTPUT_NAME=simple_repeats
#     /global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv}

#     for treatment in H3K27ac H3K36me3 H3K27me3; do
#         TREATMENT_BAM=/global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/sorted.${treatment}_1.sam.bam
#         CONTROL_BAM=/global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/sorted.Input_1.sam.bam
#         /global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $TREATMENT_BAM -d $CONTROL_BAM -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${treatment} -s ${sv}
#     done

#     DENSITY_BED=gc
#     SUB_OUTPUT_NAME=gc
#     /global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv}

#     DENSITY_BED=guy11_fungap_out_12_28_20.justgenes.bed
#     SUB_OUTPUT_NAME=genes
#     /global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv}

#     ## could probably add SCOs here or effectors or different types of genes as well
# done

sv=DEL.no_region_scale

OUTPUT_NAME=plot_profile_all_${sv}_no_survivor

REGIONS_BED=all_${sv}.bed
GENOME_FILE=guy11_genome_baoetal2017_with_70-15_mito.fasta
THREADS=${SLURM_NTASKS}
WINDOWS=10

source activate /global/scratch/users/pierrj/conda_envs/deeptools

DENSITY_BED=guy11_genome_baoetal2017_with_70-15_mito.fasta.out.gff
SUB_OUTPUT_NAME=tes
/global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv}

DENSITY_BED=guy11_genome_baoetal2017_with_70-15_mito.fasta.simplerepeats.bed
SUB_OUTPUT_NAME=simple_repeats
/global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv}

for treatment in H3K27ac H3K36me3 H3K27me3; do
    TREATMENT_BAM=/global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/sorted.${treatment}_1.sam.bam
    CONTROL_BAM=/global/scratch/users/pierrj/eccDNA/pipeline_tests/histone_marks/samples_merged/sorted.Input_1.sam.bam
    /global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $TREATMENT_BAM -d $CONTROL_BAM -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${treatment} -s ${sv}
done

DENSITY_BED=gc
SUB_OUTPUT_NAME=gc
/global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv}

DENSITY_BED=guy11_fungap_out_12_28_20.justgenes.bed
SUB_OUTPUT_NAME=genes
/global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv}