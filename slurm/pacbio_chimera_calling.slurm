#!/bin/bash
#SBATCH --job-name=pacbio_chimera_calling
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


# cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio/

## ccs demux pacbio sequences in all three SMRT cells then sort by barcode
# barcodes="/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/barcodes.fasta"
# while read sample; do sbatch -A ac_kvkallow --job-name=$sample.demux_ccs --export=sample=$sample,barcodes=$barcodes /global/home/users/pierrj/git/slurm/css_demux_mass_submission.slurm ; done < mapfile_smartcells

# while read barcode; do
#     echo ${barcode}
#     bamtools merge -in SMRTcell_1/SMRTcell_1_demux.AJI_${barcode}--OF_PacF${barcode}.bam -in SMRTcell_2/SMRTcell_2_demux.AJI_${barcode}--OF_PacF${barcode}.bam -in SMRTcell_3/SMRTcell_3_demux.AJI_${barcode}--OF_PacF${barcode}.bam -out SMRTcell_all/merged_demux.AJI_${barcode}--OF_PacF${barcode}.bam
# done < barcodelist

# mv merged_demux.AJI_17--OF_PacF17.bam G3_1A_pacbio.sam
# mv merged_demux.AJI_18--OF_PacF18.bam G3_1C_pacbio.sam
# mv merged_demux.AJI_19--OF_PacF19.bam G3_2A_pacbio.sam
# mv merged_demux.AJI_20--OF_PacF20.bam G3_2B_pacbio.sam
# mv merged_demux.AJI_21--OF_PacF21.bam G3_2C_pacbio.sam
# mv merged_demux.AJI_22--OF_PacF22.bam G3_3A_pacbio.sam
# mv merged_demux.AJI_9--OF_PacF9.bam G3_1B_pacbio.sam
# mv merged_demux.AJI_11--OF_PacF11.bam G3_3B_pacbio.sam

# 

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio/SMRTcell_all

while read sample; do
    cp ${sample}/${sample}_pacbio.bam /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/structural_variation/chimeras
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/structural_variation/chimeras

GENE_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.renamed.bed

while read sample; do
    samtools fastq /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio/SMRTcell_all/${sample}/${sample}_pacbio.bam > ${sample}_pacbio.fastq
    /global/scratch/users/pierrj/minimap2/minimap2 -ax asm20 /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_calling_pacbio/guy11_genome_baoetal2017.mmi ${sample}_pacbio.fastq > ${sample}_aligned.bam
    /global/home/users/pierrj/git/bash/get_chimeric_eccs.sh -s ${sample} -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/uniq.filtered.sorted.${sample}.bam -p ${sample}_aligned.bam -g ${GENE_BEDFILE}
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio/SMRTcell_all/mapfile