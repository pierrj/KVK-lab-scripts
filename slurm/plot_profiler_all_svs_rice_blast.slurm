#!/bin/bash
#SBATCH --job-name=plot_profiler_all_svs_rice_blast
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


cd /global/scratch/users/pierrj/PAV_SV/SV/dels_re_rice_blast

FASTA_FILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta

## get svs
if [ -f all_svs_list ]; then
    rm all_svs_list
fi
while read SAMPLE; do
    echo /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/${SAMPLE}/all/${SAMPLE}.all.vcf >> all_svs_list
done < /global/scratch/users/pierrj/sv_calling_moryzae/run_1_24_2022_guy11/gladieux_accessions

## merge
# in order #
# distance between breakpoints
# number of supporting callers
# take type into account
# take strand into account
# disabled?
# minimum size
# output
/global/scratch/users/pierrj/SURVIVOR/Debug/SURVIVOR merge all_svs_list 1000 1 1 1 0 0 all_svs.vcf
/global/scratch/users/pierrj/SURVIVOR/Debug/SURVIVOR vcftobed all_svs.vcf 1 10000000000000000 all_svs.bed

## remove mitochondria
mv all_svs.bed all_svs.bed.old
grep -v 'Supercontig_7.9' all_svs.bed.old > all_svs.bed


## get TE annotation file
if [ ! -f "guy11_genome_baoetal2017.fasta.out.gff" ]; then
    module load repeatmasker/4.1.1
    RM_LIB=/global/scratch/users/pierrj/PAV_SV/PAV/te_annotations_anne/anne_repbase_te_lib_fixed
    RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . $FASTA_FILE
fi

## get the profile plots
source activate /global/scratch/users/pierrj/conda_envs/deeptools

sv=DEL
OUTPUT_NAME=plot_profile_all_DEL
REGIONS_BED=all_DEL.bed
GENOME_FILE=guy11_genome_baoetal2017.fasta
THREADS=${SLURM_NTASKS}
WINDOWS=10

awk -v sv=$sv -v OFS='\t' '{ if ( $11 == sv ) {print $12, $13, $15}}' all_svs.bed > all_${sv}.bed

DENSITY_BED=guy11_genome_baoetal2017.fasta.out.gff
SUB_OUTPUT_NAME=tes
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

SUB_OUTPUT_NAME=cm_expression
MAPFILE_FOR_AVERAGING=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11/zhang_et_al_2019_complete_medium_bamfile_mapfile

/global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d `cat $MAPFILE_FOR_AVERAGING | tr '\n' ' ' | xargs | sed -e 's/ / -d /g'` \
    -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv}

SUB_OUTPUT_NAME=ip_expression
MAPFILE_FOR_AVERAGING=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11/zhang_et_al_2019_in_planta_bamfile_mapfile

/global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d `cat $MAPFILE_FOR_AVERAGING | tr '\n' ' ' | xargs | sed -e 's/ / -d /g'` \
    -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv}

SUB_OUTPUT_NAME=eccdna
MAPFILE_FOR_AVERAGING=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11_eccdnas/bamfile_mapfile

/global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d `cat $MAPFILE_FOR_AVERAGING | tr '\n' ' ' | xargs | sed -e 's/ / -d /g'` \
    -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv}

## to get mapfiles for averageing ##
# while read SRA; do realpath ${SRA}.Aligned.sortedByCoord.out.bam; done < $SRA_LIST > ${SRA_LIST}_bamfile_mapfile
# while read sample; do realpath ${sample}.bam; done < $MAPFILE > bamfile_mapfile