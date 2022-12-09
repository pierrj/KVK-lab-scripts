#!/bin/bash
#SBATCH --job-name=plot_profiler_all_svs_wheat_blast
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

WORKING_DIR=/global/scratch/users/pierrj/PAV_SV/SV/dels_wheat_blast

cd $WORKING_DIR

FASTA_FILE=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_to_annotate/GCA_004785725.2_ASM478572v2_genomic.fna
FASTA_FILE_TRIMMED=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_to_annotate/GCA_004785725.2_ASM478572v2_genomic.no_mito_no_mini.fna

module load seqtk

grep '>' $FASTA_FILE | grep chromosome | grep -v mini | awk '{print substr($0,2)}' | seqtk subseq ${FASTA_FILE} - > $FASTA_FILE_TRIMMED

FASTA_FILE=$FASTA_FILE_TRIMMED

basename_fasta_file=$(basename $FASTA_FILE)

SV_DIR=/global/scratch/users/pierrj/sv_calling_moryzae/wheat_blast_11_2_22
SV_MAPFILE=wheat_blast_illumina_seq_runs.txt

## get svs
if [ -f all_svs_list ]; then
    rm all_svs_list
fi
while read SAMPLE; do
    echo ${SV_DIR}/${SAMPLE}/all/${SAMPLE}.all.vcf >> all_svs_list
done < ${SV_DIR}/${SV_MAPFILE}

## merge
# in order #
# distance between breakpoints
# number of supporting callers
# take type into account
# take strand into account
# disabled?
# minimum size
# output
/global/scratch/users/pierrj/SURVIVOR/Debug/SURVIVOR merge all_svs_list 100 1 1 1 0 0 all_svs.vcf
/global/scratch/users/pierrj/SURVIVOR/Debug/SURVIVOR vcftobed all_svs.vcf 1 10000000000000000 all_svs.bed

## remove mitochondria and mini
mv all_svs.bed all_svs.bed.old
grep -v 'CP060337.1' all_svs.bed.old | grep -v 'CP060338.1' > all_svs.bed


## get TE annotation file
if [ ! -f "${basename_fasta_file}.out.gff" ]; then
    module load repeatmasker/4.1.1
    RM_LIB=/global/scratch/users/pierrj/PAV_SV/PAV/te_annotations_anne/anne_repbase_te_lib_fixed
    RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . $FASTA_FILE
fi

## get the profile plots
source activate /global/scratch/users/pierrj/conda_envs/deeptools

## generate some necessary files
samtools faidx ${FASTA_FILE}
cut -f1-2 ${FASTA_FILE}.fai > ${FASTA_FILE}.chromsizes

sv=DEL
GENOME_FILE=$FASTA_FILE
GENOME_CHROMSIZES=${FASTA_FILE}.chromsizes
THREADS=${SLURM_NTASKS}
WINDOWS=10
GENE_BED=/global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all/random_forest/all_beds_fixed/GCA_004785725.2_ASM478572v2.genes.bed
TE_GFF=${basename_fasta_file}.out.gff

grep -v 'CP060337.1' $GENE_BED | grep -v 'CP060338.1' > b71_gene_bed_no_mini_no_mito.bed

GENE_BED=b71_gene_bed_no_mini_no_mito.bed

awk -v sv=$sv -v OFS='\t' '{ if ( $11 == sv ) {print $12, $13, $15}}' all_svs.bed > all_${sv}.bed
bedtools intersect -a all_${sv}.bed -b $GENE_BED -wa -u -f 0.5 > gene_${sv}.bed
bedtools intersect -a all_${sv}.bed -b $TE_GFF -wa -u -f 0.5 > te_${sv}.bed

for i in $(seq 1 100); do cat all_${sv}.bed; done > all_${sv}.bed.replicated

bedtools shuffle -i all_${sv}.bed.replicated -g $GENOME_CHROMSIZES > shuffled_${sv}.bed

awk -v OFS='\t' '{print $1, $4, $5}' $TE_GFF | grep -v "#" > all_tes.bed

cp $GENE_BED all_genes.bed

for REGIONS_BED in all_${sv}.bed gene_${sv}.bed te_${sv}.bed shuffled_${sv}.bed all_genes.bed all_tes.bed; do
    OUTPUT_NAME=plot_profile_${REGIONS_BED}

    DENSITY_BED=$TE_GFF
    SUB_OUTPUT_NAME=tes
    FLANKING_DIST=5000
    /global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv} -f $FLANKING_DIST

    DENSITY_BED=$GENE_BED
    SUB_OUTPUT_NAME=genes
    FLANKING_DIST=5000
    /global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o ${OUTPUT_NAME}_${SUB_OUTPUT_NAME} -s ${sv} -f $FLANKING_DIST
done