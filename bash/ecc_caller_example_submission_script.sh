#!/bin/bash
mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11bao_7015mito_spikeins_tweaked_osativa_worganelles_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm; done < mapfile



mapfile="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.contignames"
genome_bwa="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GRCh37.p13_bwa"
sample="SRR6315407"
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm

mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_w_organelles_bwa"
sample="OsCT1"
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm

mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
sample="OsSR"
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_nomap.slurm

mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_w_organelles_bwa"
sample="OsCT1"
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm

mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_w_organelles_bwa"
sample="OsCT1"
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm

mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm; done < mapfile


mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
sample=ERR2660591
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_inplace.slurm


mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_with_spikeins_tweaked_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm; done < mapfile_CT_PQ


mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
while read sample; do sbatch --job-name=$sample.rarefaction --export=sample=$sample,mapfile=$mapfile /global/home/users/pierrj/git/slurm/rarefaction_plot_mass_submission.slurm; done < mapfile


mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11bao_7015mito_spikeins_tweaked_osativa_worganelles_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm; done < mapfile_IF

sample=IF_1A
mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11bao_7015mito_spikeins_tweaked_osativa_worganelles_bwa"
sbatch --job-name=$sample.ecc_caller -p savio --ntasks-per-node=20 --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm


mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_w_organelles_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm; done < mapfile

mapfile="/global/scratch/users/pierrj/references/ORSA_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_w_organelles_bwa"
sample=RC_2B
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm


mapfile="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm; done < mapfile

mapfile="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_bwa"
sample=RC_2C
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm

mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
sample=G3_1A
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm

mapfile="/global/scratch/users/pierrj/references/moryzae_70-15_ref_with_mito.contignames"
genome_bwa="/global/scratch/users/pierrj/references/moryzae_70-15_ref_with_mito_bwa"
sample=G3_1A
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_inplace.slurm

mapfile="/global/scratch/users/pierrj/references/moryzae_70-15_ref_with_mito.contignames"
genome_bwa="/global/scratch/users/pierrj/references/moryzae_70-15_ref_with_mito_bwa"
sample=ERR2660591
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_inplace.slurm

mapfile="/global/scratch/users/pierrj/ragoo_runs/G3/70-15_nobreak/70-15_nobreak.ragoo.contignames"
genome_bwa="/global/scratch/users/pierrj/ragoo_runs/G3/70-15_nobreak/70-15_nobreak.ragoo_bwa"
sample=G3_1A_70-15_ragoo
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm

mapfile="/global/scratch/users/pierrj/ragoo_runs/G3/70-15/70-15.ragoo.contignames"
genome_bwa="/global/scratch/users/pierrj/ragoo_runs/G3/70-15/70-15.ragoo_bwa"
sample=G3_1A_70-15_ragoo_breaks
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome.slurm

mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
sample=G3_1A_nomapq0
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_nomapq0.slurm


mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_nomap.slurm; done < mapfile

mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_withmapq0.slurm; done < mapfile

mapfile="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_bwa"
while read sample; do sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_withmapq0.slurm; done < mapfile

mapfile="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.contignames"
genome_bwa="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GRCh37.p13_bwa"
sample=1000DCAT
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_withmapq0.slurm

mapfile="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.contignames"
genome_bwa="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GRCh37.p13_bwa"
sample=glp-1
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_withmapq0.slurm

mapfile="/global/scratch/users/pierrj/eccDNA/2018_moller/references/Celegans_WBcel235_genomic.contignames"
genome_bwa="/global/scratch/users/pierrj/eccDNA/2018_moller/references/Celegans_WBcel235_complete_bwa"
sample=L1GexoII
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_withmapq0.slurm


mapfile="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_bwa"
sample=RC_3A
sbatch --job-name=$sample.ecc_caller --partition=savio2_bigmem --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_withmapq0.slurm

mapfile="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_bwa"
sample=RC_3C
sbatch --job-name=$sample.ecc_caller --partition=savio2_bigmem --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_withmapq0.slurm

mapfile="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"
genome_bwa="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_bwa"
sample=ERR2660591
sbatch --job-name=$sample.ecc_caller --export=sample=$sample,mapfile=$mapfile,genome_bwa=$genome_bwa /global/home/users/pierrj/git/slurm/ecc_caller_anygenome_withmapq0.slurm