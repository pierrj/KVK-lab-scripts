#!/bin/bash
#SBATCH --job-name=specific_region_splitreads
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
cd /global/scratch/users/pierrj/genome/2018_foster
# sample=G3_1A
# awk '$1=="MQOP01000002.1" && $2 > 1996800 && $2 < 2004730 {print $4}' exactlytwice.all.mergedandpe.${sample}_bwamem.bed > targetregion_readlist
# touch targetregion_distant_splitreads
# while read splitread;
# do
#     grep ${splitread} exactlytwice.all.mergedandpe.${sample}_bwamem.bed  >> targetregion_distant_splitreads
# done < targetregion_readlist
# awk '( ($1!="MQOP01000002.1") || ($1=="MQOP01000002.1" && ($2 < 1996800 || $2 > 2004730)))' targetregion_distant_splitreads > targetregion_distant_splitreads.bed
# bedtools getfasta -fi /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta -bed targetregion_distant_splitreads.bed > targetregion_distant_splitreads.fasta
# blastn -query targetregion_distant_splitreads.fasta -db /global/scratch/users/pierrj/references/ty1_copia_scaf2_1997kb.fasta -out targetregion_distant_splitreads.blastn -num_threads 20 -outfmt "6 qseqid qlen sseqid slen pident nident score qstart qend sstart send evalue" -max_target_seqs 1

# start_LTR1=1,996,945
# end_LTR1=1,997,507
# start_LTR2=2,004,067
# end_LTR2=2,004,610

sample=G3_1A
awk '( ($1=="MQOP01000002.1" && $2 > 1996945 && $2 < 1997507) || ($1=="MQOP01000002.1" && $2 > 2004067 && $2 < 2004610) ) {print $4}' exactlytwice.all.mergedandpe.${sample}_bwamem.bed > targetregion_LTRs_readlist
touch targetregion_LTRs_distant_splitreads
while read splitread;
do
    grep ${splitread} exactlytwice.all.mergedandpe.${sample}_bwamem.bed  >> targetregion_LTRs_distant_splitreads
done < targetregion_LTRs_readlist
awk '( ($1!="MQOP01000002.1") || ($1=="MQOP01000002.1" && ($2 < 1996945 || $2 > 2004610)))' targetregion_LTRs_distant_splitreads > targetregion_LTRs_distant_splitreads.bed
bedtools getfasta -fi /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta -bed targetregion_distant_splitreads.bed > targetregion_distant_splitreads.fasta
blastn -query targetregion_distant_splitreads.fasta -db /global/scratch/users/pierrj/references/ty1_copia_scaf2_1997kb.fasta -out targetregion_distant_splitreads.blastn -num_threads 20 -outfmt "6 qseqid qlen sseqid slen pident nident score qstart qend sstart send evalue" -max_target_seqs 1
