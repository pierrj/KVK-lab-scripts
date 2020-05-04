#!/bin/bash
#SBATCH --job-name=te_comparisons_bygroup
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/tecomparisons
table=groupcount.libnormalized
while read sample;
do
    awk '{ seen[$9] += $17} END { for (i in seen) print seen[i] }' ${sample}.teintersect.readsperkb.libnormalized.bed > ${sample}.${table}
done < /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/mapfile
while read sample;
do
    awk '{ seen[$9] += $17} END { for (i in seen) print seen[i] }' ${sample}.teintersect.readsperkb.libnormalized.bed > ${sample}.${table}
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
sample=CT_1A
awk '{ seen[$9] += $17} END { for (i in seen) print i }' ${sample}.${table} | awk -F":" '{print (NF>1)? $NF : ""}' > ${table}.firstcolumn
ls -1 *.${table} > ${table}.list
paste ${table}.firstcolumn $(cat list | tr "\n" " ") > ${table}.table.old
ls -1 *.${table} | awk '{print substr($0, 0, 5)}' > ${table}.first_row.old
cat <(echo x) ${table}.first_row.old | tr "\n" "\t" | sed '$ s/.$//' > ${table}.first_row
cat ${table}.first_row <(echo) ${table}.table.old > ${table}.table
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
}' ${table}.table > ${table}.table.transposed


table=groupcount.spikenormalized
while read sample;
do
    awk '{ seen[$9] += $17} END { for (i in seen) print seen[i] }' ${sample}.teintersect.readsperkb.libnormalized.bed > ${sample}.${table}
done < /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/mapfile
sample=CT_1A
awk '{ seen[$9] += $17} END { for (i in seen) print i }' ${sample}.${table} | awk -F":" '{print (NF>1)? $NF : ""}' > ${table}.firstcolumn
ls -1 *.${table} > ${table}.list
paste ${table}.firstcolumn $(cat list | tr "\n" " ") > ${table}.table.old
ls -1 *.${table} | awk '{print substr($0, 0, 5)}' > ${table}.first_row.old
cat <(echo x) ${table}.first_row.old | tr "\n" "\t" | sed '$ s/.$//' > ${table}.first_row
cat ${table}.first_row <(echo) ${table}.table.old > ${table}.table
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
}' ${table}.table > ${table}.table.transposed