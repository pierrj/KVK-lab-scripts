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
CHROM_SIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes
LTR_TE_LOCATIONS=/global/scratch/pierrj/references/te_annotations/moryzae/moryzae_copiaandgypsy_locs.bed

# cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot/large_ecc_density_old .
# cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot/micro_dna_density_old .

# awk -v OFS='\t' '{print $1, $2-100, $2, $3}' large_ecc_density_old > large_ecc_density_old.bg
# awk -v OFS='\t' '{print $1, $2-100, $2, $3}' micro_dna_density_old > micro_dna_density_old.bg

# bedtools merge -i large_ecc_density_old.bg -c 4 -o mean -d -1 > large_ecc_density_old_new.bg
# bedGraphToBigWig large_ecc_density_old_new.bg ${GENOME_CHROMSIZES} large_ecc_density_old.bw

# computeMatrix scale-regions -p ${SLURM_NTASKS} -S large_ecc_density_old.bw \
#                             -R ${GENE_GTF} \
#                             --beforeRegionStartLength 10000 \
#                             --regionBodyLength 1000 \
#                             --afterRegionStartLength 10000 \
#                             --skipZeros -o large_eccdna.mat.gz

# plotProfile -m large_eccdna.mat.gz \
#             -out large_eccdna.pdf \
#             --numPlotsPerRow 1 \
#             --plotTitle "large" \
#             --outFileNameData large_eccdna.tab

# bedtools merge -i micro_dna_density_old.bg -c 4 -o mean -d -1 > micro_dna_density_old_new.bg
# bedGraphToBigWig micro_dna_density_old_new.bg ${GENOME_CHROMSIZES} micro_dna_density_old.bw

# computeMatrix scale-regions -p ${SLURM_NTASKS} -S micro_dna_density_old.bw \
#                             -R ${GENE_GTF} \
#                             --beforeRegionStartLength 10000 \
#                             --regionBodyLength 1000 \
#                             --afterRegionStartLength 10000 \
#                             --skipZeros -o micro_dna.mat.gz

# plotProfile -m micro_dna.mat.gz \
#             -out micro_dna.pdf \
#             --numPlotsPerRow 1 \
#             --plotTitle "micro" \
#             --outFileNameData micro_dna.tab


## simulated




track=micro_dna_density
group=micro_dnas

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample
do
cd ${sample}
bedtools shuffle -i ${group}_splitreads -g ${GENOME_CHROMSIZES} -excl ${LTR_TE_LOCATIONS} > ${group}_splitreads.permuted
bedtools coverage -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot/guy11_genome_baoetal2017.100bpwindows \
    -b ${group}_splitreads.permuted -g ${CHROM_SIZES} | \
    awk -v OFS='\t' '{print $1, $3, $4}' > ${sample}.${group}.permuted.genomecoverage.bed
cd ..
done < mapfile

if [ -f "${group}.normalization_table" ]; then
    rm ${group}.normalization_table
fi

while read sample; do
    cd ${sample}
    total=$(wc -l ${group}_splitreads.permuted | awk '{print $1/1000000}')
    echo -e ${sample}'\t'${total} >> ../${group}.normalization_table
    cd ..
done < mapfile

/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -m sample_mapfile \
    -t ${sample}.${group}.permuted.genomecoverage.bed -n ${group}.normalization_table -y t

/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1 -b 1 -c 3 -n n

awk -v OFS='\t' '{print $1, $2, $3}' G3.normalized_binned > /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot/${track}_old_permuted

cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/circos_plot/micro_dna_density_old_permuted .

awk -v OFS='\t' '{print $1, $2-100, $2, $3}' micro_dna_density_old_permuted > micro_dna_density_old_permuted.bg


bedtools merge -i micro_dna_density_old_permuted.bg -c 4 -o mean -d -1 > micro_dna_density_old_permuted.bg
bedGraphToBigWig micro_dna_density_old_permuted.bg ${GENOME_CHROMSIZES} micro_dna_density_old_permuted.bw

computeMatrix scale-regions -p ${SLURM_NTASKS} -S micro_dna_density_old_permuted.bw \
                            -R ${GENE_GTF} \
                            --beforeRegionStartLength 1000 \
                            --regionBodyLength 1000 \
                            --afterRegionStartLength 1000 \
                            --skipZeros -o micro_dna_density_old_permuted.mat.gz

plotProfile -m micro_dna_density_old_permuted.mat.gz \
            -out micro_dna_density_old_permuted.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "micro" \
            --outFileNameData micro_dna_density_old_permuted.tab