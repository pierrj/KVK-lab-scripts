#!/bin/bash
#SBATCH --job-name=star_map
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
for (( i = 47 ; i < 50; i++))
do
    cd /global/scratch/users/pierrj/rnaseq/moryzae/2012_soanes/${i}_star
    STAR --runThreadN 24 \
        --genomeDir /global/scratch/users/pierrj/references/GUY11_PacBio_merge_KROJ_renamed_star \
        --readFilesIn /global/scratch/users/pierrj/rnaseq/moryzae/2012_soanes/SRR2988${i}.sra.fastq \
        --quantMode GeneCounts
    awk 'NR>4 {print $3}' ReadsPerGene.out.tab > ../${i}.ReadsPerGene.out.tab
    cd ..
    normalize_factor=$(awk '{SUM+=$1}END{print SUM/1000000}' ${i}.ReadsPerGene.out.tab)
    awk -v N=$normalize_factor '{print $1/N}' ${i}.ReadsPerGene.out.tab > normalized.${i}.ReadsPerGene.out.tab
done

cd /global/scratch/users/pierrj/rnaseq/moryzae/2012_soanes/
awk 'NR>4 {print $1}' 47_star/ReadsPerGene.out.tab > genecount_firstcolumn
paste genecount_firstcolumn normalized.47.ReadsPerGene.out.tab normalized.48.ReadsPerGene.out.tab normalized.49.ReadsPerGene.out.tab > genecount_table
awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' genecount_table > genecount_table_average.old
paste genecount_firstcolumn genecount_table_average.old > genecount_table_average