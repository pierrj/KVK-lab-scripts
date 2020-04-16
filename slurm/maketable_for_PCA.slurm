#!/bin/bash
#SBATCH --job-name=maketable_for_wholegenome_PCA
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina
while read sample;
do
    cd ${sample}
    total=$(samtools view -c -F 2048 sorted.mergedandpe.${sample}_bwamem.bam)
    echo -e $sample'\t'$total >> /global/scratch/users/pierrj/eccDNA/pipeline_tests/pca/illumina_libsizes
    cp genomecoverage.mergedandpe.${sample}_bwamem.bed /global/scratch/users/pierrj/eccDNA/pipeline_tests/pca/
    cd ..
done < mapfile
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress
while read sample;
do
    cd ${sample}
    total=$(samtools view -c -F 2048 sorted.mergedandpe.${sample}_bwamem.bam)
    echo -e $sample'\t'$total >> /global/scratch/users/pierrj/eccDNA/pipeline_tests/pca/illumina_libsizes
    cp genomecoverage.mergedandpe.${sample}_bwamem.bed /global/scratch/users/pierrj/eccDNA/pipeline_tests/pca/
    cd ..
done < mapfile
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/pca
ls genomecoverage.mergedandpe.*_bwamem.bed  | awk '{print substr($0,28,5)}' | sort | uniq > mapfile
while read sample; 
do
    normalize_factor=$(grep ${sample} illumina_libsizes | awk '{print $2/1000000}')
    awk -v N=$normalize_factor '{sum+=$3} NR%100==0 {print sum/100/N; sum =0}' genomecoverage.mergedandpe.${sample}_bwamem.bed > binned_libnormalized.genomecoverage.mergedandpe.${sample}_bwamem.bed
done < mapfile
normalize_factor=$(grep CT_1A illumina_libsizes | awk '{print $2/1000000}')
awk -v OFS='-' -v N=$normalize_factor '{sum+=$3} NR%100==0 {print $1,$2; sum =0}' genomecoverage.mergedandpe.CT_1A_bwamem.bed > first_column
ls -1 binned_libnormalized.genomecoverage.mergedandpe.${sample}_bwamem.bed > list
paste first_column $(cat list | tr "\n" " ") > binned_libnormalized.coverage.mergedandpe.table.old
ls -1 binned_libnormalized.genomecoverage.mergedandpe.*_bwamem.bed | awk '{print substr($0, 49, 5)}' > first_row.old
cat <(echo x) first_row.old | tr "\n" "\t" | sed '$ s/.$//' > first_row
cat first_row <(echo) binned_libnormalized.coverage.mergedandpe.table.old > binned_libnormalized.coverage.mergedandpe.table
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
}' binned_libnormalized.coverage.mergedandpe.table > transposed_table.txt