#!/bin/bash
#SBATCH --job-name=gather_data_for_rf_prediction
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


## actual variables
PROJ_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/test_rf_data_generation/mg8
GENOME_OUTPUT_STRING=MG8
PROTEIN_FASTA=${PROJ_DIR}/Magnaporthe_oryzae.MG8.pep.all.fa
PROTEIN_GFF=${PROJ_DIR}/Magnaporthe_oryzae.MG8.55.gff3
GENOME_FASTA=${PROJ_DIR}/Magnaporthe_oryzae.MG8.dna.toplevel.fa

## shouldnt change much
PFAM_DIR=/global/scratch/users/pierrj/pfam/
RM_LIB=/global/scratch/users/pierrj/PAV_SV/PAV/te_annotations_anne/anne_repbase_te_lib_fixed
GUY11_PROTEOME=/global/scratch/users/pierrj/PAV_SV/PAV/test_rf_data_generation/guy11_data/GUY11_fungap_out_prot_filtered.faa
GUY11_DATA_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/test_rf_data_generation/guy11_data

cd $PROJ_DIR

## some prep
awk -v OFS='\t' '{if ($3 == "gene") {print $1, $4, $5, $9}}' $PROTEIN_GFF > ${GENOME_OUTPUT_STRING}_gene.bed
PROTEIN_BED=${PROJ_DIR}/${GENOME_OUTPUT_STRING}_gene.bed
samtools faidx ${GENOME_FASTA}
cut -f1-2 ${GENOME_FASTA}.fai > ${GENOME_FASTA}.chromsizes
GENOME_CHROM_SIZES=${GENOME_FASTA}.chromsizes

### GATHER DATA
if [ -f slurm_ids ]; then
    rm slurm_ids
fi
## get pfam annotations
sbatch --parsable \
--export=ALL,PROJ_DIR=$PROJ_DIR,GENOME_OUTPUT_STRING=$GENOME_OUTPUT_STRING,PROTEIN_FASTA=$PROTEIN_FASTA,PROTEIN_GFF=$PROTEIN_GFF,GENOME_FASTA=$GENOME_FASTA,PFAM_DIR=$PFAM_DIR,RM_LIB=$RM_LIB,GUY11_PROTEOME=$GUY11_PROTEOME,GUY11_DATA_DIR=$GUY11_DATA_DIR,PROTEIN_BED=$PROTEIN_BED,GENOME_CHROM_SIZES=$GENOME_CHROM_SIZES \
                    /global/home/users/pierrj/git/slurm/pfam_scan_per_proteome.slurm >> slurm_ids

## get GC
sbatch --parsable -A co_minium \
--export=ALL,PROJ_DIR=$PROJ_DIR,GENOME_OUTPUT_STRING=$GENOME_OUTPUT_STRING,PROTEIN_FASTA=$PROTEIN_FASTA,PROTEIN_GFF=$PROTEIN_GFF,GENOME_FASTA=$GENOME_FASTA,PFAM_DIR=$PFAM_DIR,RM_LIB=$RM_LIB,GUY11_PROTEOME=$GUY11_PROTEOME,GUY11_DATA_DIR=$GUY11_DATA_DIR,PROTEIN_BED=$PROTEIN_BED,GENOME_CHROM_SIZES=$GENOME_CHROM_SIZES \
/global/home/users/pierrj/git/slurm/gc_gene_stats.slurm >> slurm_ids

## repeatmask
sbatch --parsable \
--export=ALL,PROJ_DIR=$PROJ_DIR,GENOME_OUTPUT_STRING=$GENOME_OUTPUT_STRING,PROTEIN_FASTA=$PROTEIN_FASTA,PROTEIN_GFF=$PROTEIN_GFF,GENOME_FASTA=$GENOME_FASTA,PFAM_DIR=$PFAM_DIR,RM_LIB=$RM_LIB,GUY11_PROTEOME=$GUY11_PROTEOME,GUY11_DATA_DIR=$GUY11_DATA_DIR,PROTEIN_BED=$PROTEIN_BED,GENOME_CHROM_SIZES=$GENOME_CHROM_SIZES \
/global/home/users/pierrj/git/slurm/repeatmask_per_genome.slurm >> slurm_ids

## effector domains
sbatch --parsable -A co_minium \
--export=ALL,PROJ_DIR=$PROJ_DIR,GENOME_OUTPUT_STRING=$GENOME_OUTPUT_STRING,PROTEIN_FASTA=$PROTEIN_FASTA,PROTEIN_GFF=$PROTEIN_GFF,GENOME_FASTA=$GENOME_FASTA,PFAM_DIR=$PFAM_DIR,RM_LIB=$RM_LIB,GUY11_PROTEOME=$GUY11_PROTEOME,GUY11_DATA_DIR=$GUY11_DATA_DIR,PROTEIN_BED=$PROTEIN_BED,GENOME_CHROM_SIZES=$GENOME_CHROM_SIZES \
 /global/home/users/pierrj/git/slurm/effector_predictor.slurm >> slurm_ids

## GO
sbatch --parsable -A co_minium \
--export=ALL,PROJ_DIR=$PROJ_DIR,GENOME_OUTPUT_STRING=$GENOME_OUTPUT_STRING,PROTEIN_FASTA=$PROTEIN_FASTA,PROTEIN_GFF=$PROTEIN_GFF,GENOME_FASTA=$GENOME_FASTA,PFAM_DIR=$PFAM_DIR,RM_LIB=$RM_LIB,GUY11_PROTEOME=$GUY11_PROTEOME,GUY11_DATA_DIR=$GUY11_DATA_DIR,PROTEIN_BED=$PROTEIN_BED,GENOME_CHROM_SIZES=$GENOME_CHROM_SIZES \
/global/home/users/pierrj/git/slurm/pannzer_per_proteome.slurm >> slurm_ids

## xfer info from guy11
sbatch --parsable -A co_minium \
--export=ALL,PROJ_DIR=$PROJ_DIR,GENOME_OUTPUT_STRING=$GENOME_OUTPUT_STRING,PROTEIN_FASTA=$PROTEIN_FASTA,PROTEIN_GFF=$PROTEIN_GFF,GENOME_FASTA=$GENOME_FASTA,PFAM_DIR=$PFAM_DIR,RM_LIB=$RM_LIB,GUY11_PROTEOME=$GUY11_PROTEOME,GUY11_DATA_DIR=$GUY11_DATA_DIR,PROTEIN_BED=$PROTEIN_BED,GENOME_CHROM_SIZES=$GENOME_CHROM_SIZES \
/global/home/users/pierrj/git/slurm/transfer_info_from_guy11_per_genome.slurm >> slurm_ids

## make sure everything finishes
if [ -f slurm_ids ]; then
    # wait for all jobs to finish
    sbatch -W --dependency=afterok:$(cat slurm_ids | head -c -1 | tr '\n' ':') /global/home/users/pierrj/git/slurm/dummy_job.slurm

    wait
fi

## bring it all together now
