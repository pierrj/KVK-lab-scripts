#!/bin/bash
#SBATCH --job-name=blast_eccDNA_ends_100bp
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/homology/

if [ -f "homology_table_10000_shuffled" ]; then
    rm homology_table_10000_shuffled
fi

while read eccdna;
do
    echo $eccdna > eccdna
    ## problem here with eccs at start and ends of scaffolds but oh well
    ## problem here with really small eccs as well
    ## problem here with TE eccDNAs
    awk -v OFS='\t' '{ print $1, $2-10, $2+10}' eccdna > eccdna_start
    awk -v OFS='\t' '{ print $1, $3-10, $3+10}' eccdna > eccdna_end
    bedtools getfasta -fi /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta -fo eccdna_start.fasta -bed eccdna_start
    bedtools getfasta -fi /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta -fo eccdna_end.fasta -bed eccdna_end
    blastn -query eccdna_start.fasta -subject eccdna_end.fasta -word_size 4 -task blastn-short -max_hsps 1 -perc_identity 100 -outfmt  "6 qseqid sseqid qstart qend sstart send length" > blastoutput ## add -word_size 4 here and test
    if [ -f blastoutput ]
    then
        echo -e $(awk '{print $7}' blastoutput ) >> homology_table_10000_shuffled
        rm blastoutput
    else
        echo -e "0" >> homology_table_10000_shuffled
    fi
done < 10000_random.bed

## blastn -query eccDNA_start.fasta -subject eccDNA_end.fasta -task blastn-short -perc_identity 100