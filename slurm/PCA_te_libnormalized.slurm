#!/bin/bash
#SBATCH --job-name=PCA_te_libnormalized
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
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina
while read sample;
do
    cd ${sample}
    samtools view -b sorted.mergedandpe.${sample}_bwamem.bam $(cat /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames | tr "\n" " ") > guy11.sorted.allmapq.mergedandpe.${sample}_bwamem.bam
    bedtools intersect -sorted -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes -F 1 -wa -c -a /global/scratch/users/pierrj/references/GUY11_pacbio_ET_KROJ.features.forintersect.gff -b guy11.sorted.allmapq.mergedandpe.${sample}_bwamem.bam > /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.bed
    awk -v OFS='\t' '{print $1, $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$15/$13*1000}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.bed > /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.readsperkb.bed
    normalize_factor=$(awk '{sum += $16} END {print sum/100000}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.readsperkb.bed)
    awk -v N=$normalize_factor -v OFS='\t' '{print $1, $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$16/N}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.readsperkb.bed > /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.readsperkb.libnormalized.bed
    awk '{print $17}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.readsperkb.libnormalized.bed > /global/scratch/users/pierrj/eccDNA/pipeline_tests/pca/te_libnormalized/${sample}.teintersect.forpca.bed
    cd ..
done < mapfile
cd /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress
while read sample;
do
    cd ${sample}
    samtools view -b sorted.mergedandpe.${sample}_bwamem.bam $(cat /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames | tr "\n" " ") > guy11.sorted.allmapq.mergedandpe.${sample}_bwamem.bam
    bedtools intersect -sorted -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes -F 1 -wa -c -a /global/scratch/users/pierrj/references/GUY11_pacbio_ET_KROJ.features.forintersect.gff -b guy11.sorted.allmapq.mergedandpe.${sample}_bwamem.bam > /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.bed
    awk -v OFS='\t' '{print $1, $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$15/$13*1000}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.bed > /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.readsperkb.bed
    normalize_factor=$(awk '{sum += $16} END {print sum/100000}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.readsperkb.bed)
    awk -v N=$normalize_factor -v OFS='\t' '{print $1, $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$16/N}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.readsperkb.bed > /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.readsperkb.libnormalized.bed
    awk '{print $17}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.readsperkb.libnormalized.bed > /global/scratch/users/pierrj/eccDNA/pipeline_tests/pca/te_libnormalized/${sample}.teintersect.forpca.bed
    cd ..
done < mapfile
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/pca/te_libnormalized
awk -v OFS=';' '{print $1,$4,$5,$9}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/CT_1A.teintersect.readsperkb.libnormalized.bed > first_column
ls -1 *.teintersect.forpca.bed > list
paste first_column $(cat list | tr "\n" " ") > libnormalized.te.table.old
ls -1 *.teintersect.forpca.bed | awk '{print substr($0, 0, 5)}' > first_row.old
cat <(echo x) first_row.old | tr "\n" "\t" | sed '$ s/.$//' > first_row
cat first_row <(echo) libnormalized.te.table.old > libnormalized.te.table
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
}' libnormalized.te.table > libnormalized.te.table.transposed
ls *.teintersect.forpca.bed  | awk '{print substr($0,0,4)}' | sort | uniq > techrepsmerged/mapfile
while read sample; 
do
    paste $(find . -maxdepth 1 -name "${sample}*.teintersect.forpca.bed" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > techrepsmerged/${sample}.teintersect.forpca.bed
done < techrepsmerged/mapfile
cp first_column techrepsmerged/first_column
cd techrepsmerged
ls -1 *.teintersect.forpca.bed > list
paste first_column $(cat list | tr "\n" " ") > libnormalized.te.techrepsmerged.table.old
ls -1 *.teintersect.forpca.bed | awk '{print substr($0, 0, 4)}' > first_row.old
cat <(echo x) first_row.old | tr "\n" "\t" | sed '$ s/.$//' > first_row
cat first_row <(echo) libnormalized.te.techrepsmerged.table.old > libnormalized.te.techrepsmerged.table
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
}' libnormalized.te.techrepsmerged.table > libnormalized.te.techrepsmerged.table.transposed