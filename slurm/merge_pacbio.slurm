#!/bin/bash
#SBATCH --job-name=merge_pacbio
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio

while read barcode; do
    echo ${barcode}
    bamtools merge -in SMRTcell_1/SMRTcell_1_demux.AJI_${barcode}--OF_PacF${barcode}.bam -in SMRTcell_2/SMRTcell_2_demux.AJI_${barcode}--OF_PacF${barcode}.bam -in SMRTcell_3/SMRTcell_3_demux.AJI_${barcode}--OF_PacF${barcode}.bam -out SMRTcell_all/merged_demux.AJI_${barcode}--OF_PacF${barcode}.bam
    samtools view -c SMRTcell_all/merged_demux.AJI_${barcode}--OF_PacF${barcode}.bam
done < barcodelist

SMRTcell_1_demux.AJI_18--OF_PacF19.bam