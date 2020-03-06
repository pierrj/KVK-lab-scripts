#!/bin/bash
#SBATCH --job-name=doesntmatter
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
cd $SLURM_SUBMIT_DIR
mkdir ${sample}
cd ${sample}
find ../ -name "${sample}*R1*" | xargs -I '{}' mv {} ./${sample}_R1.fastq
find ../ -name "${sample}*R2*" | xargs -I '{}' mv {} ./${sample}_R2.fastq
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o tmp.trimmed.${sample}_R1.fastq -p tmp.trimmed.${sample}_R2.fastq ${sample}_R1.fastq ${sample}_R2.fastq
/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/SeqPrep/SeqPrep/SeqPrep -f tmp.trimmed.${sample}_R1.fastq -r tmp.trimmed.${sample}_R2.fastq -1 tmp.seqprep.trimmed.${sample}_R1.fastq -2 tmp.seqprep.trimmed.${sample}_R2.fastq -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -s tmp.merged.trimmed.${sample}
bwa mem -t 24 /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/guy11_genome_baoetal2017_bwa tmp.seqprep.trimmed.${sample}_R1.fastq tmp.seqprep.trimmed.${sample}_R2.fastq -o tmp.seqprep.trimmed.${sample}_bwamem.sam
bwa mem -t 24 /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/guy11_genome_baoetal2017_bwa tmp.merged.trimmed.${sample} -o tmp.merged.trimmed.${sample}_bwamem.sam
samtools view -S -b tmp.seqprep.trimmed.${sample}_bwamem.sam > tmp.seqprep.trimmed.${sample}_bwamem.bam
samtools view -S -b tmp.merged.trimmed.${sample}_bwamem.sam > tmp.merged.trimmed.${sample}_bwamem.bam
bamtools merge -in tmp.seqprep.trimmed.${sample}_bwamem.bam -in tmp.merged.trimmed.${sample}_bwamem.bam -out tmp.mergedandpe.${sample}_bwamem.bam

## get splitread input file
samtools view -f 81 -F 4 tmp.mergedandpe.${sample}_bwamem.bam > tmp.reverseread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 145 -F 4 tmp.mergedandpe.${sample}_bwamem.bam > tmp.reverseread2.mergedandpe.${sample}_bwamem.sam
samtools view -f 65 -F 20 tmp.mergedandpe.${sample}_bwamem.bam > tmp.forwardread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 129 -F 20 tmp.mergedandpe.${sample}_bwamem.bam > tmp.forwardread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.reverseread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.reverseread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.forwardread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread2.mergedandpe.${sample}_bwamem.sam tmp.forwardread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam
cat  tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam > samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam

samtools view tmp.mergedandpe.${sample}_bwamem.bam | awk '{ if (($2 == 83 || $2 == 147 ) && $9 > 0) print $0 ; else if (($2 == 99 || $2 == 163) && $9 <0) print $0}' | cat <(samtools view -H tmp.mergedandpe.${sample}_bwamem.bam) - | samtools view -b -h - > tmp.outwardfacing.mergedandpe.${sample}_bwamem.bam
bedtools bamtobed -i tmp.outwardfacing.mergedandpe.${sample}_bwamem.bam | sort -k 4,4 > outwardfacing.mergedandpe.${sample}_bwamem.bed

samtools sort tmp.mergedandpe.${sample}_bwamem.bam -o tmp.sorted.mergedandpe.${sample}_bwamem.bam
bedtools genomecov -d -ibam tmp.sorted.mergedandpe.${sample}_bwamem.bam > genomecoverage.mergedandpe.${sample}_bwamem.bed

rm tmp.*

ipcluster start -n $SLURM_NTASKS --cluster-id="slurm-${SLURM_JOBID}" &
sleep 45
ipython /global/home/users/pierrj/git/python/ecc_caller.py samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam outwardfacing.mergedandpe.${sample}_bwamem.bed genomecoverage.mergedandpe.${sample}_bwamem.bed ${sample}
ipcluster stop --cluster-id="slurm-${SLURM_JOBID}"