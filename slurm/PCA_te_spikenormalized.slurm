#!/bin/bash
#SBATCH --job-name=PCA_te_spikenormalized
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
    normalize_factor=$(grep ${sample} /global/scratch/users/pierrj/eccDNA/pipeline_tests/normalization/SRcount_slopes | awk '{print $2}')
    awk -v N=$normalize_factor -v OFS='\t' '{print $1, $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$16/N}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.readsperkb.bed > /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.readsperkb.spikenormalized.bed
    awk '{print $17}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/${sample}.teintersect.readsperkb.spikenormalized.bed > /global/scratch/users/pierrj/eccDNA/pipeline_tests/pca/te_spikenormalized/${sample}.teintersect.forpca.bed
    cd ..
done < mapfile
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/pca/te_spikenormalized
awk -v OFS=';' '{print $1,$4,$5,$9}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons/CT_1A.teintersect.readsperkb.spikenormalized.bed > first_column
ls -1 *.teintersect.forpca.bed > list
paste first_column $(cat list | tr "\n" " ") > spikenormalized.te.table.old
ls -1 *.teintersect.forpca.bed | awk '{print substr($0, 0, 5)}' > first_row.old
cat <(echo x) first_row.old | tr "\n" "\t" | sed '$ s/.$//' > first_row
cat first_row <(echo) spikenormalized.te.table.old > spikenormalized.te.table
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
}' spikenormalized.te.table > spikenormalized.te.table.transposed
ls *.teintersect.forpca.bed  | awk '{print substr($0,0,4)}' | sort | uniq > techrepsmerged/mapfile
while read sample; 
do
    paste $(find . -maxdepth 1 -name "${sample}*.teintersect.forpca.bed" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > techrepsmerged/${sample}.teintersect.forpca.bed
done < techrepsmerged/mapfile
cp first_column techrepsmerged/first_column
cd techrepsmerged
ls -1 *.teintersect.forpca.bed > list
paste first_column $(cat list | tr "\n" " ") > spikenormalized.te.techrepsmerged.table.old
ls -1 *.teintersect.forpca.bed | awk '{print substr($0, 0, 4)}' > first_row.old
cat <(echo x) first_row.old | tr "\n" "\t" | sed '$ s/.$//' > first_row
cat first_row <(echo) spikenormalized.te.techrepsmerged.table.old > spikenormalized.te.techrepsmerged.table
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
}' spikenormalized.te.techrepsmerged.table > spikenormalized.te.techrepsmerged.table.transposed