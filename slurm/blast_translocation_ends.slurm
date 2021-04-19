#!/bin/bash
#SBATCH --job-name=blast_translocation_ends_100bp
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/processing_translocations

if [ -f "homology_table_100" ]; then
    rm homology_table_100
fi

while read translocation;
do
    echo $translocation > translocation
    awk -v OFS='\t' '{ print $1, $2-50, $2+50}' translocation > translocation_start
    awk -v OFS='\t' '{ print $1, $3-50, $3+50}' translocation > translocation_end
    bedtools getfasta -fi /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta -fo translocation_start.fasta -bed translocation_start
    bedtools getfasta -fi /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta -fo translocation_end.fasta -bed translocation_end
    blastn -query translocation_start.fasta -subject translocation_end.fasta -task blastn-short > blastoutput ## add -word_size 4 here and test
    cat blastoutput
    if [ -f blastoutput ]
    then
        echo -e $translocation'\t'$(awk '{print $7}' blastoutput ) >> homology_table_100
        rm blastoutput
    else
        echo -e $translocation'\t'"0" >> homology_table_100
    fi
done < sorted.translocations.bed