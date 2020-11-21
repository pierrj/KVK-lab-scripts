#!/bin/bash
cd $1
nohup python /Users/pierrj/fungap_local/FunGAP/fungap_patched.py \
  --output_dir fungap_out \
  --trans_read_1 SRR8842990_1.fastq \
  --trans_read_2 SRR8842990_2.fastq \
  --genome_assembly guy11_genome_baoetal2017.fasta  \
  --augustus_species magnaporthe_grisea  \
  --sister_proteome prot_db.faa  \
  --busco_dataset sordariomycetes_odb10 \
  --num_cores 2 >& run.out &