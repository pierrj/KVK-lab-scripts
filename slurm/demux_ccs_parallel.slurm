#!/bin/bash
#SBATCH --job-name=SMRTcell3_ccs_demux
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=24
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/SMRTcell3
cp m54213_190515_135227.subreads.bam subreads.bam
samtools view -h subreads.bam > subreads.sam
split --number=l/24 -d --additional-suffix=.sam subreads.sam linesplit
samtools view -H subreads.sam > head.linesplit0.sam
for (( i = 1 ; i < 24; i++)); do cp head.linesplit0.sam head.linesplit$i.sam ; done
for (( i = 0 ; i < 10; i++)); do cat head.linesplit$i.sam linesplit0$i.sam > split$i.subreads.sam ; done
for (( i = 10 ; i < 24; i++)); do cat head.linesplit$i.sam linesplit$i.sam > split$i.subreads.sam ; done
for (( i = 0 ; i < 24; i++)); do samtools view -b split$i.subreads.sam > split$i.subreads.bam; done 
source activate python=2.7
for (( i = 0 ; i < 24; i++))
do
    srun --ntasks=1 --nodes=1 ccs --reportFile=ccs_report_split$i.txt split$i.subreads.bam $i.split.ccs.bam &
done
wait
sleep 60
bamtools merge -in 0.split.ccs.bam -in 1.split.ccs.bam -in 2.split.ccs.bam -in 3.split.ccs.bam -in 4.split.ccs.bam -in 5.split.ccs.bam -in 6.split.ccs.bam -in 7.split.ccs.bam -in 8.split.ccs.bam -in 9.split.ccs.bam -in 10.split.ccs.bam -in 11.split.ccs.bam -in 12.split.ccs.bam -in 13.split.ccs.bam -in 14.split.ccs.bam -in 15.split.ccs.bam -in 16.split.ccs.bam -in 17.split.ccs.bam -in 18.split.ccs.bam -in 19.split.ccs.bam -in 20.split.ccs.bam -in 21.split.ccs.bam -in 22.split.ccs.bam -in 23.split.ccs.bam -out ccs.bam
lima --different --dump-removed --split-bam-named --ccs ccs.bam /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/barcodes.fasta demux.bam