#!/bin/bash
#SBATCH --job-name=WGS-allSRs_no_orientation
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/genome/2018_foster/
sample=ERR2660591
file="mergedandpe.${sample}_bwamem.bam"

samtools view -f 65 -F 4 ${file} > mergedandpe.1.${sample}_bwamem.sam
samtools view -f 129 -F 4 ${file} > mergedandpe.2.${sample}_bwamem.sam

awk 'NR==FNR{a[$1]++; next} a[$1]==2' mergedandpe.1.${sample}_bwamem.sam mergedandpe.1.${sample}_bwamem.sam > exactlytwice.mergedandpe.1.${sample}_bwamem.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' mergedandpe.2.${sample}_bwamem.sam mergedandpe.2.${sample}_bwamem.sam > exactlytwice.mergedandpe.2.${sample}_bwamem.sam
cat  exactlytwice.mergedandpe.1.${sample}_bwamem.sam exactlytwice.mergedandpe.2.${sample}_bwamem.sam > exactlytwice.all.mergedandpe.${sample}_bwamem.bam

qualityfilter_file="exactlytwice.all.mergedandpe.${sample}_bwamem.bam"


awk '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if((a !~ /[DMIHS]/ && int(a) > 19 ) || (b !~ /[DMIHS]/ && int(b) > 19)) print $0}' ${qualityfilter_file} > qualityfiltered.${qualityfilter_file}
awk 'NR==FNR{a[$1]++; next} a[$1]==2' qualityfiltered.${qualityfilter_file} qualityfiltered.${qualityfilter_file} > exactlytwice.qualityfiltered.${qualityfilter_file}
samtools view -b -h <(cat <(samtools view -H /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/mergedandpe.${sample}_bwamem.bam) exactlytwice.qualityfiltered.${qualityfilter_file}) > SRs.no_orientation.mergedandpe.${sample}_bwamem.bam
bedtools bamtobed -i SRs.no_orientation.mergedandpe.${sample}_bwamem.bam | sort -k4,4 -k2,2n > SRs.no_orientation.mergedandpe.${sample}_bwamem.bed

samtools sort SRs.no_orientation.mergedandpe.${sample}_bwamem.bam > sorted.SRs.no_orientation.mergedandpe.${sample}_bwamem.bam
samtools index sorted.SRs.no_orientation.mergedandpe.${sample}_bwamem.bam