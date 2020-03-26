#!/bin/bash
#SBATCH --job-name=getSRs_pBR322
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
# cd $SLURM_SUBMIT_DIR
# mkdir ${sample}
# cd ${sample}
# find ../ -name "${sample}*R1*" | xargs -I '{}' mv {} ./${sample}_R1.fastq
# find ../ -name "${sample}*R2*" | xargs -I '{}' mv {} ./${sample}_R2.fastq
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/normalization/guy11_and_spikeins
sample=CT_1A_pBR322
## get splitread input file
samtools view -f 81 -F 4 mergedandpe.${sample}_bwamem.bam > tmp.reverseread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 145 -F 4 mergedandpe.${sample}_bwamem.bam > tmp.reverseread2.mergedandpe.${sample}_bwamem.sam
samtools view -f 65 -F 20 mergedandpe.${sample}_bwamem.bam > tmp.forwardread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 129 -F 20 mergedandpe.${sample}_bwamem.bam > tmp.forwardread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.reverseread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.reverseread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.forwardread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread2.mergedandpe.${sample}_bwamem.sam tmp.forwardread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam
cat  tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam > samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam

rm tmp.*

awk '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if((a !~ /[DMIHS]/ && int(a) > 19 ) || (b !~ /[DMIHS]/ && int(b) > 19)) print $0}' samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam > qualityfiltered.samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam
samtools view -b -h <(cat <(samtools view -H mergedandpe.${sample}_bwamem.bam) exactlytwice.qualityfiltered.samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam) > SRs.mergedandpe.${sample}_bwamem.bam
bedtools bamtobed -i SRs.mergedandpe.${sample}_bwamem.bam | sort -k 4,4 -k 2,2 > SRs.mergedandpe.${sample}_bwamem.bed