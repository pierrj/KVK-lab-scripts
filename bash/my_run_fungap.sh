#!/bin/bash
cd $1
genome=$(basename $1)
echo start >> /Users/pierrj/fungap_runs/moryzae_genomes/annotation_run_tracker/${genome}_run_tracker
date +"%T" >> /Users/pierrj/fungap_runs/moryzae_genomes/annotation_run_tracker/${genome}_run_tracker

python /Users/pierrj/fungap_local/FunGAP/fungap_patched.py \
  --output_dir fungap_out \
  --trans_read_1 SRR8842990_1.fastq \
  --trans_read_2 SRR8842990_2.fastq \
  --genome_assembly ${genome}_genomic.fna  \
  --augustus_species magnaporthe_grisea  \
  --sister_proteome prot_db.faa  \
  --busco_dataset sordariomycetes_odb10 \
  --num_cores 1 >& run.out

if [ ! -f ${1}/fungap_out/fungap_out/fungap_out.gff3 ]; then
    echo $1 >> /Users/pierrj/fungap_runs/moryzae_genomes/annotation_run_tracker/did_not_complete
fi

echo end >> /Users/pierrj/fungap_runs/moryzae_genomes/annotation_run_tracker/${genome}_run_tracker
date +"%T" >> /Users/pierrj/fungap_runs/moryzae_genomes/annotation_run_tracker/${genome}_run_tracker