#!/bin/bash
#SBATCH --job-name=PCA_gene_spikenormalized
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
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress
while read sample;
do
    cd ${sample}
    bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/eccDNA/pipeline_tests/sp_enrichment/moryzae_all_trimmed.uniquefiltered.gff -b /global/scratch/users/pierrj/eccDNA/pipeline_tests/merging/parallel.confirmed.${sample}.bed | awk -v OFS='\t' '{print $4, $5}' > /global/scratch/users/pierrj/eccDNA/pipeline_tests/pca/gene_spikenormalized/genes.parallel.confirmed.${sample}.bed
    normalize_factor=$(grep ${sample} /global/scratch/users/pierrj/eccDNA/pipeline_tests/normalization/SRcount_slopes | awk '{print $2}')
    awk -v N=$normalize_factor '{print $2/N; sum =0}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/pca/gene_spikenormalized/genes.parallel.confirmed.${sample}.bed > /global/scratch/users/pierrj/eccDNA/pipeline_tests/pca/gene_spikenormalized/spikenormalized.genes.parallel.confirmed.${sample}.bed
    cd ..
done < mapfile
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/pca/gene_spikenormalized
awk '{print $1}' spikenormalized.genes.parallel.confirmed.CT_1A.bed > first_column
ls -1 spikenormalized.genes.parallel.confirmed.*.bed > list
paste first_column $(cat list | tr "\n" " ") > spikenormalized.gene.table.old
ls -1 spikenormalized.genes.parallel.confirmed.*.bed | awk '{print substr($0, 41, 5)}' > first_row.old
cat <(echo x) first_row.old | tr "\n" "\t" | sed '$ s/.$//' > first_row
cat first_row <(echo) spikenormalized.gene.table.old > spikenormalized.gene.table
awk 'BEGIN { FS=OFS="\t" }
{
    for (rowNr=1;rowNr<=NF;rowNr++) {
        cell[rowNr,NR] = $rowNr
    }
    maxRows = (NF > maxRows ? NF : maxRows)
    maxCols = NR
}
END {
    for (rowNr=1;rowNr<=maxRows;rowNr++) {
        for (colNr=1;colNr<=maxCols;colNr++) {
            printf "%s%s", cell[rowNr,colNr], (colNr < maxCols ? OFS : ORS)
        }
    }
}' spikenormalized.gene.table > spikenormalized.gene.table.transposed
mkdir techrepsmerged
ls spikenormalized.genes.parallel.confirmed.*.bed  | awk '{print substr($0,41,4)}' | sort | uniq > techrepsmerged/mapfile
while read sample; 
do
    paste $(find . -maxdepth 1 -name "spikenormalized.genes.parallel.confirmed.${sample}*" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > techrepsmerged/spikenormalized.genes.parallel.confirmed.${sample}.bed
done < techrepsmerged/mapfile
cp first_column techrepsmerged/first_column
cd techrepsmerged
ls -1 spikenormalized.genes.parallel.confirmed.*.bed > list
paste first_column $(cat list | tr "\n" " ") > spikenormalized.gene.techrepsmerged.table.old
ls -1 spikenormalized.genes.parallel.confirmed.*.bed | awk '{print substr($0, 41, 4)}' > first_row.old
cat <(echo x) first_row.old | tr "\n" "\t" | sed '$ s/.$//' > first_row
cat first_row <(echo) spikenormalized.gene.techrepsmerged.table.old > spikenormalized.gene.techrepsmerged.table
awk 'BEGIN { FS=OFS="\t" }
{
    for (rowNr=1;rowNr<=NF;rowNr++) {
        cell[rowNr,NR] = $rowNr
    }
    maxRows = (NF > maxRows ? NF : maxRows)
    maxCols = NR
}
END {
    for (rowNr=1;rowNr<=maxRows;rowNr++) {
        for (colNr=1;colNr<=maxCols;colNr++) {
            printf "%s%s", cell[rowNr,colNr], (colNr < maxCols ? OFS : ORS)
        }
    }
}' spikenormalized.gene.techrepsmerged.table > spikenormalized.gene.techrepsmerged.table.transposed