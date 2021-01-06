#!/bin/bash
#SBATCH --job-name=blast_eccDNA_ends_100bp
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/homology/
rm homology_table_100
touch homology_table_100
while read eccDNA;
do
    echo $eccDNA > eccDNA
    awk -v OFS='\t' '{ print $1, $2-50, $2+50, $4}' eccDNA > eccDNA_start
    awk -v OFS='\t' '{ print $1, $3-50, $3+50, $4}' eccDNA > eccDNA_end
    bedtools getfasta -fi /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta -fo eccDNA_start.fasta -bed eccDNA_start
    bedtools getfasta -fi /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta -fo eccDNA_end.fasta -bed eccDNA_end
    blastn -query eccDNA_start.fasta -subject eccDNA_end.fasta -task blastn-short -strand 'plus' -max_hsps 1 -perc_identity 100 -outfmt  "6 qseqid sseqid qstart qend sstart send length" > blastoutput ## add -word_size 4 here and test
    if [ -s blastoutput ]
    then
        echo -e $eccDNA'\t'$(awk '{print $7}' blastoutput ) >> homology_table_100
    else
        echo -e $eccDNA'\t'"0" >> homology_table_100
    fi
    rm blastoutput
done < ecccaller_output.G3_1A.trimmed.bed

## blastn -query eccDNA_start.fasta -subject eccDNA_end.fasta -task blastn-short -perc_identity 100
