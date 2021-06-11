


echo /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina_if_v_g3/ >> mapfile_ecc_count_dir

echo /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta >> mapfile_ecc_count_genome_file

echo /global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.gff3 >> mapfile_ecc_count_gff_file

echo /global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_copia_locs.bed >> mapfile_ecc_count_copiafile

echo /global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_gypsy_locs.bed >> mapfile_ecc_count_gypsyfile

echo mapfile >> mapfile_ecc_count_mapfile

paste mapfile_ecc_count_dir mapfile_ecc_count_genome_file mapfile_ecc_count_gff_file mapfile_ecc_count_copiafile mapfile_ecc_count_gypsyfile mapfile_ecc_count_mapfile > mapfile_ecc_count