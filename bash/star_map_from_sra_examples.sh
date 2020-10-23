#!/bin/bash
ADAPTER_ONE=AGATCGGAAGAGCACACGTCTGAAC
ADAPTER_TWO=AGATCGGAAGAGCGTCGTGTAGGGA
sample=SRR7642295
GENOME_DB=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_star
sbatch --job-name=$sample.starmap --export=sample=$sample,GENOME_DB=$GENOME_DB,ADAPTER_ONE=$ADAPTER_ONE,ADAPTER_TWO=$ADAPTER_TWO /global/home/users/pierrj/git/slurm/star_map_from_sra_pe.slurm


sample=SRR7642295
GENOME_DB=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_star
sbatch --job-name=$sample.starmap --export=sample=$sample,GENOME_DB=$GENOME_DB /global/home/users/pierrj/git/slurm/star_map_from_sra_pe.slurm

sample=SRR8842990
GENOME_DB=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_star
sbatch --job-name=$sample.starmap_pe --export=sample=$sample,GENOME_DB=$GENOME_DB /global/home/users/pierrj/git/slurm/star_map_from_sra_pe.slurm

sample=SRR306516
GENOME_DB=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_star
sbatch --job-name=$sample.starmap_se --export=sample=$sample,GENOME_DB=$GENOME_DB /global/home/users/pierrj/git/slurm/star_map_from_sra_se.slurm

