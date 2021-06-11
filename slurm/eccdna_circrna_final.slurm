#!/bin/bash
#SBATCH --job-name=eccdna_circrna
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/eccdna_gene_expression/zhang_et_al_data/zhang_et_al_rnaseq_metadata.txt .
cp /global/scratch/pierrj/eccDNA/pipeline_tests/eccdna_gene_expression/zhang_et_al_data/zhang_et_al_rnaseq_accessions.txt .

GENOME_FILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta

THREADS=${SLURM_NTASKS}

while read SRA; do
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SRA} -O .
    /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SRA}.sra
done < zhang_et_al_rnaseq_accessions.txt

awk -v FS=',' -v OFS='\t' '{print $1,$50}' zhang_et_al_rnaseq_metadata.txt > metadata_biorep.txt
awk -v OFS='\t' '{print $2, $3}' metadata_biorep.txt | tail -n+2 | sort | uniq > bioreps
while read biorep ; do
awk -v b=$biorep '{if ($2==b) {print $1".sra.fastq"}}' metadata_biorep.txt > accessions
cat $(cat accessions | tr "\n" " ") > zhang_${biorep}.fastq
done < bioreps

while read biorep; do
    bwa mem -T 19 -t ${SLURM_NTASKS} ${GENOME_FILE} zhang_${biorep}.fastq -o zhang_${biorep}.sam
    perl /global/scratch/users/pierrj/ciri2/CIRI_v2.0.5/CIRI_v2.0.5.pl -T ${SLURM_NTASKS} -I zhang_${biorep}.sam -O zhang_${biorep}.ciri_out -F ${GENOME_FILE}
done < bioreps

cat zhang_1.ciri_out zhang_2.ciri_out zhang_3.ciri_out | awk -v OFS='\t' '{print $2, $3, $4}' | sort -k1,1 -k2,2n | uniq | head -n-1 > zhang_all.ciri_out

if [ -f "all.ecc_caller_out.details.nolowq.txt" ]; then
    rm all.ecc_caller_out.details.nolowq.txt
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/


while read sample; do
cd ${sample}
    cat ${sample}.ecc_caller_out.details.nolowq.txt >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas/all.ecc_caller_out.details.nolowq.txt
cd ..
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

sort -k1,1 -k2,2n all.ecc_caller_out.details.nolowq.txt | awk -v OFS='\t' '{print $1, $2, $3}' | uniq > noduplicates.all.ecc_caller_out.details.nolowq.txt