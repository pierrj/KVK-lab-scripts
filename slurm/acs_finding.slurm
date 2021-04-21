#!/bin/bash
#SBATCH --job-name=acs_finding
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

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/acs_finding

module load meme

if [ -f permute_fimo_out.txt ]; then
    rm permute_fimo_out.txt
fi
for i in {1..100}; do 
    echo $i
    if [ -d shuffled_fimo ]; then
        rm -r shuffled_fimo
    fi
    bedtools shuffle -i notes.sorted.uniq.all.confirmedsplitreads.bed -g guy11_genome_baoetal2017.chromsizes -excl guy11_genome_baoetal2017.fasta.out.gff > shuffled
    bedtools getfasta -fi guy11_genome_baoetal2017.fasta -bed shuffled > shuffled.fasta
    fimo --verbosity 1 -o shuffled_fimo acs_meme.txt shuffled.fasta
    awk '$6>=17' shuffled_fimo/fimo.txt | awk '{print $2}' | sort | uniq | wc -l >> permute_fimo_out.txt
done