

while read sample
do
mkdir ${sample}
cd ${sample}
cp /global/scratch/pierrj/eccDNA/magnaporthe_pureculture/illumina_w_merge/${sample}/${sample}_R1.fastq .
cp /global/scratch/pierrj/eccDNA/magnaporthe_pureculture/illumina_w_merge/${sample}/${sample}_R2.fastq .
cd ..
done < mapfile


while read sample
do
cd ${sample}
rm *.fastq
rm *.bam
rm *.sam
cd ..
done < mapfile


while read sample
do
cd ${sample}
cp ecccaller_output.${sample}.renamed.bed /global/scratch/users/pierrj/eccDNA/2018_moller/full_run_prenomerge/
cp ecccaller_output.${sample}.renamed.details.tsv /global/scratch/users/pierrj/eccDNA/2018_moller/full_run_prenomerge/
cd ..
done < mapfile

while read sample
do
cd ${sample}
cp ecccaller_output.${sample}.renamed.bed /global/scratch/users/pierrj/eccDNA/2020_wang/outputs_pre_nomerge
cp ecccaller_output.${sample}.renamed.details.tsv /global/scratch/users/pierrj/eccDNA/2020_wang/outputs_pre_nomerge
cd ..
done < mapfile

for i in $(seq 8752175 1 8752206)
do
    scancel $i
done

while read sample
do
cd ${sample}
cp ecccaller_output.${sample}.renamed.bed /global/scratch/users/pierrj/eccDNA/2017_lanciano/full_run/outputs_pre_nomerge
cp ecccaller_output.${sample}.renamed.details.tsv /global/scratch/users/pierrj/eccDNA/2017_lanciano/full_run/outputs_pre_nomerge
cd ..
done < mapfile

/global/scratch/users/pierrj/eccDNA/2015_moller/full_run/outputs_pre_nomerge

while read sample
do
cd ${sample}
cp ecccaller_output.${sample}.renamed.bed /global/scratch/users/pierrj/eccDNA/2015_moller/full_run/outputs_pre_nomerge
cp ecccaller_output.${sample}.renamed.details.tsv /global/scratch/users/pierrj/eccDNA/2015_moller/full_run/outputs_pre_nomerge
cd ..
done < mapfile
