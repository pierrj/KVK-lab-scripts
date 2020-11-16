#!/bin/bash
echo /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/ >> mapfile_ecc_count_dir
echo /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/ >> mapfile_ecc_count_dir
echo /global/scratch/users/pierrj/eccDNA/2018_moller/full_run >> mapfile_ecc_count_dir
echo /global/scratch/users/pierrj/eccDNA/2018_moller/full_run >> mapfile_ecc_count_dir
echo /global/scratch/users/pierrj/eccDNA/2017_lanciano/full_run >> mapfile_ecc_count_dir
echo /global/scratch/users/pierrj/eccDNA/2017_lanciano/full_run >> mapfile_ecc_count_dir
echo /global/scratch/users/pierrj/eccDNA/2017_lanciano/full_run >> mapfile_ecc_count_dir
echo /global/scratch/users/pierrj/eccDNA/2017_lanciano/full_run >> mapfile_ecc_count_dir
echo /global/scratch/users/pierrj/eccDNA/2017_lanciano/full_run >> mapfile_ecc_count_dir
echo /global/scratch/users/pierrj/eccDNA/2015_moller/full_run >> mapfile_ecc_count_dir
echo /global/scratch/users/pierrj/eccDNA/2015_moller/full_run >> mapfile_ecc_count_dir
echo /global/scratch/users/pierrj/eccDNA/2015_moller/full_run >> mapfile_ecc_count_dir
echo /global/scratch/users/pierrj/eccDNA/2015_moller/full_run >> mapfile_ecc_count_dir


echo moryzae >> mapfile_ecc_count_output
echo myrice >> mapfile_ecc_count_output
echo human_muscle >> mapfile_ecc_count_output
echo human_leukocytes >> mapfile_ecc_count_output
echo ara_wt >> mapfile_ecc_count_output
echo ara_epi >> mapfile_ecc_count_output
echo ory_callus >> mapfile_ecc_count_output
echo ory_leaf >> mapfile_ecc_count_output
echo ory_seed >> mapfile_ecc_count_output
echo yeast_d >> mapfile_ecc_count_output
echo yeast_z >> mapfile_ecc_count_output
echo yeast_g >> mapfile_ecc_count_output
echo yeast_s >> mapfile_ecc_count_output

echo /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta >> mapfile_ecc_count_genome_file
echo /global/scratch/users/pierrj/references/ORSA_IRGSP-1.0.fasta >> mapfile_ecc_count_genome_file
echo /global/scratch/users/pierrj/references/GRCh37.p13.fasta >> mapfile_ecc_count_genome_file
echo /global/scratch/users/pierrj/references/GRCh37.p13.fasta >> mapfile_ecc_count_genome_file
echo /global/scratch/users/pierrj/references/TAIR10.fasta >> mapfile_ecc_count_genome_file
echo /global/scratch/users/pierrj/references/TAIR10.fasta >> mapfile_ecc_count_genome_file
echo /global/scratch/users/pierrj/references/ORSA_IRGSP-1.0.fasta >> mapfile_ecc_count_genome_file
echo /global/scratch/users/pierrj/references/ORSA_IRGSP-1.0.fasta >> mapfile_ecc_count_genome_file
echo /global/scratch/users/pierrj/references/ORSA_IRGSP-1.0.fasta >> mapfile_ecc_count_genome_file
echo /global/scratch/users/pierrj/references/Scer_S288C.fasta >> mapfile_ecc_count_genome_file
echo /global/scratch/users/pierrj/references/Scer_S288C.fasta >> mapfile_ecc_count_genome_file
echo /global/scratch/users/pierrj/references/Scer_S288C.fasta >> mapfile_ecc_count_genome_file
echo /global/scratch/users/pierrj/references/Scer_S288C.fasta >> mapfile_ecc_count_genome_file

echo /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta >> mapfile_ecc_count_gff_file
echo /global/scratch/users/pierrj/references/ORSA_IRGSP-1.0.gff >> mapfile_ecc_count_gff_file
echo /global/scratch/users/pierrj/references/GRCh37.p13.gff >> mapfile_ecc_count_gff_file
echo /global/scratch/users/pierrj/references/GRCh37.p13.gff >> mapfile_ecc_count_gff_file
echo /global/scratch/users/pierrj/references/TAIR10_GFF3_genes.gff >> mapfile_ecc_count_gff_file
echo /global/scratch/users/pierrj/references/TAIR10_GFF3_genes.gff >> mapfile_ecc_count_gff_file
echo /global/scratch/users/pierrj/references/ORSA_IRGSP-1.0.gff >> mapfile_ecc_count_gff_file
echo /global/scratch/users/pierrj/references/ORSA_IRGSP-1.0.gff >> mapfile_ecc_count_gff_file
echo /global/scratch/users/pierrj/references/ORSA_IRGSP-1.0.gff >> mapfile_ecc_count_gff_file
echo /global/scratch/users/pierrj/references/Scer_S288C.gff >> mapfile_ecc_count_gff_file
echo /global/scratch/users/pierrj/references/Scer_S288C.gff >> mapfile_ecc_count_gff_file
echo /global/scratch/users/pierrj/references/Scer_S288C.gff >> mapfile_ecc_count_gff_file
echo /global/scratch/users/pierrj/references/Scer_S288C.gff >> mapfile_ecc_count_gff_file

echo mapfile >> mapfile_ecc_count_mapfile
echo mapfile >> mapfile_ecc_count_mapfile
echo mapfile_muscle >> mapfile_ecc_count_mapfile
echo mapfile_leukocytes >> mapfile_ecc_count_mapfile
echo mapfile_ENA_ara_wt >> mapfile_ecc_count_mapfile
echo mapfile_ENA_ara_epi >> mapfile_ecc_count_mapfile
echo mapfile_ENA_ory_callus >> mapfile_ecc_count_mapfile
echo mapfile_ENA_ory_leaf >> mapfile_ecc_count_mapfile
echo mapfile_ENA_ory_seed >> mapfile_ecc_count_mapfile
echo mapfile_d >> mapfile_ecc_count_mapfile
echo mapfile_z >> mapfile_ecc_count_mapfile
echo mapfile_g >> mapfile_ecc_count_mapfile
echo mapfile_s >> mapfile_ecc_count_mapfile

paste mapfile_ecc_count_dir mapfile_ecc_count_output mapfile_ecc_count_genome_file mapfile_ecc_count_gff_file mapfile_ecc_count_mapfile > mapfile_ecc_count