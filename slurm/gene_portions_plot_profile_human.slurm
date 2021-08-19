#!/bin/bash
#SBATCH --job-name=gene_portions_plot_profile_human
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/human_plotprofile

GENOME_CHROMSIZES=/global/scratch/pierrj/references/GRCh37.p13.chromsizes
GENE_GTF=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/human/GRCh37.p13.gffread.gtf
CHROM_SIZES=/global/scratch/pierrj/references/GRCh37.p13.chromsizes
LTR_TE_LOCATIONS=/global/scratch/pierrj/references/te_annotations/human/blank
copia_file=/global/scratch/pierrj/references/te_annotations/human/blank
gypsy_file=/global/scratch/pierrj/references/te_annotations/human/blank

source activate deeptools

if [ -f "mapfile_for_normalize_and_average_filecolumn" ]; then
    rm mapfile_for_normalize_and_average_filecolumn
fi


#bedtools makewindows -g ${GENOME_CHROMSIZES} -w 10 > GRCh37.p13.10bpwindows

track=micro_dna_density
group=micro_dnas

cd /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/

while read sample
do
cd ${sample}
    bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b ${copia_file} | \
        awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
        awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $1, $2, $3}}' | \
        bedtools intersect -wao -a - -b ${gypsy_file} | \
        awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
        awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $1, $2, $3}}'| sort -k1,1 -k2,2n | awk '$3-$2 <= 400' > micro_dnas
    sed 's/[[:space:]]*$//' ${sample}.ecc_caller_out.splitreads.bed | awk 'NR==FNR{a[$0]++;next}a[$0]' micro_dnas - > micro_dnas_splitreads
    cp micro_dnas_splitreads ${sample}.micro_dnas_splitreads
    num_srs=$(wc -l ${sample}.micro_dnas_splitreads | awk '{print $1/100000}')
    bedtools coverage -a /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/human_plotprofile/GRCh37.p13.10bpwindows \
        -b ${sample}.micro_dnas_splitreads | \
        awk -v OFS='\t' -v N=${num_srs} '{print $1, $2, $3, $4/N}' > /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/human_plotprofile/${sample}.${group}.genomecoverage.bed
    echo ${sample}.${group}.genomecoverage.bed >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/human_plotprofile/mapfile_for_normalize_and_average_filecolumn
cd ..
done < mapfile_muscle

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/human_plotprofile

eval "paste $( echo $(while read sample; do
            echo "<( awk '{print \$4}' ${sample})"
        done < mapfile_for_normalize_and_average_filecolumn | tr "\n" " ") )" | \
        awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' | \
        paste <( awk '{print $1, $2, $3}' $(head -1 mapfile_for_normalize_and_average_filecolumn) ) - > ${track}.averaged

awk -v OFS='\t' '{print $1, $2, $3, $4}' ${track}.averaged > ${track}.bg

sort -k1,1 -k2,2n micro_dna_density.bg > micro_dna_density.sorted.bg

bedGraphToBigWig micro_dna_density.sorted.bg ${GENOME_CHROMSIZES} ${track}.bw

computeMatrix scale-regions -p ${SLURM_NTASKS} -S ${track}.bw \
                            -R ${GENE_GTF} \
                            --beforeRegionStartLength 2000 \
                            --regionBodyLength 1000 \
                            --afterRegionStartLength 2000 \
                            --skipZeros -o ${track}.mat.gz

plotProfile -m ${track}.mat.gz \
            -out ${track}.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "${track}" \
            --outFileNameData ${track}.tab