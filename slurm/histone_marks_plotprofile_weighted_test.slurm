#!/bin/bash
#SBATCH --job-name=histone_marks_plotprofile_weighted_test
#SBATCH --partition=savio2_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/histone_marks/weighted

if [ -f "large_eccdnas_all" ]; then
    rm large_eccdnas_all
fi

while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/large_eccdnas_splitreads >> large_eccdnas_all; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
mv large_eccdnas_all large_eccdnas_weighted

if [ -f "micro_dnas_all" ]; then
    rm micro_dnas_all
fi

while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/micro_dnas_splitreads >> micro_dnas_all; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
mv micro_dnas_all micro_dnas_weighted

cp ../large_eccdnas .
cp ../micro_dnas .
cp ../genome_baseline .

source activate deeptools

# computeMatrix scale-regions -p ${SLURM_NTASKS} -S guy11_genome_baoetal2017.50bpwindows.bw \
#                             -R large_eccdnas_weighted micro_dnas_weighted genome_baseline \
#                             --beforeRegionStartLength 1000 \
#                             --regionBodyLength 1000 \
#                             --afterRegionStartLength 1000 \
#                             --skipZeros -o gc.mat.gz

# plotProfile -m gc.mat.gz \
#             -out gc.pdf \
#             --numPlotsPerRow 1 \
#             --plotTitle "GC eccDNAs" \
#             --outFileNameData gc.tab



# computeMatrix scale-regions -p ${SLURM_NTASKS} -S guy11_genome_baoetal2017.50bpwindows.bw \
#                             -R large_eccdnas micro_dnas genome_baseline \
#                             --beforeRegionStartLength 1000 \
#                             --regionBodyLength 1000 \
#                             --afterRegionStartLength 1000 \
#                             --skipZeros -o gc_unweighted.mat.gz

# plotProfile -m gc_unweighted.mat.gz \
#             -out gc_unweighted.pdf \
#             --numPlotsPerRow 1 \
#             --plotTitle "GC eccDNAs" \
#             --outFileNameData gc_unweighted.tab



treatment=H3K27ac
bamCompare -p ${SLURM_NTASKS} -b1 ${data_dir}/sorted.${treatment}_${biorep}.sam.bam -b2 ${data_dir}/sorted.Input_${biorep}.sam.bam -o log2ratio_${treatment}_${biorep}.bw -of bigwig --scaleFactorsMethod readCount
computeMatrix scale-regions -p ${SLURM_NTASKS} -S log2ratio_${treatment}_${biorep}.bw \
                            -R large_eccdnas micro_dnas large_eccdnas_weighted micro_dnas_weighted genome_baseline \
                            --beforeRegionStartLength 300 \
                            --regionBodyLength 1000 \
                            --afterRegionStartLength 300 \
                            --skipZeros -o ${treatment}_${biorep}.mat.gz 
plotProfile -m ${treatment}_${biorep}.mat.gz \
            -out ${treatment}_${biorep}.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "${treatment} ${biorep} eccDNAs" \
            --outFileNameData ${treatment}_${biorep}.tab

treatment=H3K27me3
bamCompare -p ${SLURM_NTASKS} -b1 ${data_dir}/sorted.${treatment}_${biorep}.sam.bam -b2 ${data_dir}/sorted.Input_${biorep}.sam.bam -o log2ratio_${treatment}_${biorep}.bw -of bigwig --scaleFactorsMethod readCount
computeMatrix scale-regions -p ${SLURM_NTASKS} -S log2ratio_${treatment}_${biorep}.bw \
                            -R  large_eccdnas micro_dnas large_eccdnas_weighted micro_dnas_weighted genome_baseline \
                            --beforeRegionStartLength 300 \
                            --regionBodyLength 1000 \
                            --afterRegionStartLength 300 \
                            --skipZeros -o ${treatment}_${biorep}.mat.gz 
plotProfile -m ${treatment}_${biorep}.mat.gz \
            -out ${treatment}_${biorep}.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "${treatment} ${biorep} eccDNAs" \
            --outFileNameData ${treatment}_${biorep}.tab

treatment=H3K36me3
bamCompare -p ${SLURM_NTASKS} -b1 ${data_dir}/sorted.${treatment}_${biorep}.sam.bam -b2 ${data_dir}/sorted.Input_${biorep}.sam.bam -o log2ratio_${treatment}_${biorep}.bw -of bigwig --scaleFactorsMethod readCount
computeMatrix scale-regions -p ${SLURM_NTASKS} -S log2ratio_${treatment}_${biorep}.bw \
                            -R  large_eccdnas micro_dnas large_eccdnas_weighted micro_dnas_weighted genome_baseline \
                            --beforeRegionStartLength 300 \
                            --regionBodyLength 1000 \
                            --afterRegionStartLength 300 \
                            --skipZeros -o ${treatment}_${biorep}.mat.gz 
plotProfile -m ${treatment}_${biorep}.mat.gz \
            -out ${treatment}_${biorep}.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "${treatment} ${biorep} eccDNAs" \
            --outFileNameData ${treatment}_${biorep}.tab