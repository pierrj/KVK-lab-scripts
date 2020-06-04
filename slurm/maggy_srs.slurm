#!/bin/bash
#SBATCH --job-name=maggy_srs.slurm
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/maggy_srs/
sample=G3_1A
file="/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/guy11.repeats.sorted.mergedandpe.G3_1A_bwamem.bam"
samtools view -f 81 -F 4 ${file} > tmp.reverseread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 145 -F 4 ${file} > tmp.reverseread2.mergedandpe.${sample}_bwamem.sam
samtools view -f 65 -F 20 ${file} > tmp.forwardread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 129 -F 20 ${file} > tmp.forwardread2.mergedandpe.${sample}_bwamem.sam

##these say sam chromosome bu they aren't actually same chromosome
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.reverseread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.reverseread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.forwardread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.forwardread2.mergedandpe.${sample}_bwamem.sam tmp.forwardread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam
cat  tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam > exactlytwice.all.mergedandpe.${sample}_bwamem.bam

qualityfilter_file="exactlytwice.all.mergedandpe.${sample}_bwamem.bam"

awk '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if((a !~ /[DMIHS]/ && int(a) > 19 ) || (b !~ /[DMIHS]/ && int(b) > 19)) print $0}' ${qualityfilter_file} > qualityfiltered.${qualityfilter_file}
awk 'NR==FNR{a[$1]++; next} a[$1]==2' qualityfiltered.${qualityfilter_file} qualityfiltered.${qualityfilter_file} > exactlytwice.qualityfiltered.${qualityfilter_file}
samtools view -b -h <(cat <(samtools view -H /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/mergedandpe.${sample}_bwamem.bam) exactlytwice.qualityfiltered.${qualityfilter_file}) > SRs.mergedandpe.${sample}_bwamem.bam
bedtools bamtobed -i SRs.mergedandpe.${sample}_bwamem.bam | sort -k4,4 -k2,2n > SRs.mergedandpe.${sample}_bwamem.bed 

##to recap, these are not uniquely mapped filtered or same chromosome filtered

bedtools intersect -f 1 -wa -a SRs.mergedandpe.${sample}_bwamem.bed -b /global/scratch/users/pierrj/references/GUY11_pacbio_maggy_ltr_locs.bed > maggy.SRs.mergedandpe.${sample}_bwamem.bed
awk 'NR==FNR{a[$4]++; next} a[$4]==2' maggy.SRs.mergedandpe.${sample}_bwamem.bed maggy.SRs.mergedandpe.${sample}_bwamem.bed > exactlytwice.maggy.SRs.mergedandpe.${sample}_bwamem.bed