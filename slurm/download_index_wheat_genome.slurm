#!/bin/bash
#SBATCH --job-name=download_index_wheat_genome
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --account=co_minium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

# cd /global/scratch/projects/vector_kvklab/references/genome_fasta

# wget http://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz

# gunzip Triticum_aestivum.IWGSC.dna.toplevel.fa.gz

# cd /global/scratch/projects/vector_kvklab/references/gff

# wget http://ftp.ensemblgenomes.org/pub/plants/release-54/gff3/triticum_aestivum/Triticum_aestivum.IWGSC.54.gff3.gz

# gunzip Triticum_aestivum.IWGSC.54.gff3.gz

cd /global/scratch/projects/vector_kvklab/references/STAR_indices

mkdir Triticum_aestivum.IWGSC_plus_B71_starindex

module load STAR

WHEAT_GENOME_FASTA=/global/scratch/projects/vector_kvklab/references/genome_fasta/Triticum_aestivum.IWGSC.dna.toplevel.fa
WHEAT_GFF_FILE=/global/scratch/projects/vector_kvklab/references/gff/Triticum_aestivum.IWGSC.54.gff3

MORYZAE_GENOME_FASTA=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_to_annotate/GCA_004785725.2_ASM478572v2_genomic.fna
MORYZAE_GFF_FILE=/global/scratch/users/pierrj/fungap_runs/wheat_blast/all_gffs_fixed/GCA_004785725.2_ASM478572v2_fungap_out.fixed.gff3

cat $WHEAT_GENOME_FASTA $MORYZAE_GENOME_FASTA > /global/scratch/projects/vector_kvklab/references/genome_fasta/chinese_spring_w_b71.fasta
cat $WHEAT_GFF_FILE $MORYZAE_GFF_FILE > /global/scratch/projects/vector_kvklab/references/gff/chinese_spring_w_b71.gff3

BOTH_FASTA=/global/scratch/projects/vector_kvklab/references/genome_fasta/chinese_spring_w_b71.fasta
BOTH_GFF=/global/scratch/projects/vector_kvklab/references/gff/chinese_spring_w_b71.gff3

## index genome for STAR
STAR --runThreadN ${SLURM_NTASKS} --runMode genomeGenerate --genomeDir Triticum_aestivum.IWGSC_plus_B71_starindex \
    --genomeFastaFiles $BOTH_FASTA \
    --sjdbGTFfile $BOTH_GFF \
    --sjdbOverhang 100 \
    --genomeSAindexNbases 14 \
    --sjdbGTFtagExonParentTranscript ID \
    --sjdbGTFtagExonParentGene Parent