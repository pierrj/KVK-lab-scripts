#!/bin/bash
#SBATCH --job-name=gene_portions_plot_profile
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions

source activate deeptools

GENOME_CHROMSIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes
GENE_GTF=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.gtf

cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot/large_ecc_density_old .
cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot/micro_dna_density_old .

awk -v OFS='\t' '{print $1, $2-100, $2, $3}' large_ecc_density_old > large_ecc_density_old.bg
awk -v OFS='\t' '{print $1, $2-100, $2, $3}' micro_dna_density_old > micro_dna_density_old.bg

bedtools merge -i large_ecc_density_old.bg -c 4 -o mean -d -1 > large_ecc_density_old_new.bg
bedGraphToBigWig large_ecc_density_old_new.bg ${GENOME_CHROMSIZES} large_ecc_density_old.bw

computeMatrix scale-regions -p ${SLURM_NTASKS} -S large_ecc_density_old.bw \
                            -R ${GENE_GTF} \
                            --beforeRegionStartLength 1000 \
                            --regionBodyLength 1000 \
                            --afterRegionStartLength 1000 \
                            --skipZeros -o large_eccdna.mat.gz

plotProfile -m large_eccdna.mat.gz \
            -out large_eccdna.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "large" \
            --outFileNameData large_eccdna.tab

bedtools merge -i micro_dna_density_old.bg -c 4 -o mean -d -1 > micro_dna_density_old_new.bg
bedGraphToBigWig micro_dna_density_old_new.bg ${GENOME_CHROMSIZES} micro_dna_density_old.bw

computeMatrix scale-regions -p ${SLURM_NTASKS} -S micro_dna_density_old.bw \
                            -R ${GENE_GTF} \
                            --beforeRegionStartLength 1000 \
                            --regionBodyLength 1000 \
                            --afterRegionStartLength 1000 \
                            --skipZeros -o micro_dna.mat.gz

plotProfile -m micro_dna.mat.gz \
            -out micro_dna.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "micro" \
            --outFileNameData micro_dna.tab