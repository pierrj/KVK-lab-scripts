#!/bin/bash
#SBATCH --job-name=assemble_final_eccdna_files
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina

output_file=this_study_moryzae_illumina_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pacbio

output_file=this_study_moryzae_pacbio_eccdnas.txt

## need to fixx this

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3}' ${sample}_pacbio_eccs.bed >> ${output_file}
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/pacbio/SMRTcell_all/mapfile

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}


cd /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control

output_file=this_study_osativa_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}


cd /global/scratch/users/pierrj/eccDNA/2018_moller/full_run

output_file=2018_moller_hsapiens_muscle_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_muscle

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}


output_file=2018_moller_hsapiens_leukocyte_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_leukocytes

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}


cd /global/scratch/users/pierrj/eccDNA/2020_wang/

output_file=2021_wang_athaliana_flower_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_flower

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}

output_file=2021_wang_athaliana_leaf_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_leaf

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}

output_file=2021_wang_athaliana_root_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_root

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}

output_file=2021_wang_athaliana_stem_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_stem

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}


cd /global/scratch/users/pierrj/eccDNA/2017_lanciano/full_run/

output_file=2017_lanciano_athaliana_wt_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_ENA_ara_wt

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}

output_file=2017_lanciano_athaliana_epi12_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_ENA_ara_epi

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}

output_file=2017_lanciano_osativa_callus_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_ENA_ory_callus

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}

output_file=2017_lanciano_osativa_leaf_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_ENA_ory_leaf

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}

output_file=2017_lanciano_osativa_seed_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_ENA_ory_seed

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}

cd /global/scratch/users/pierrj/eccDNA/2015_moller/full_run/

if [ -f "sample_to_accesions_yeast" ]; then
    rm sample_to_accesions_yeast
fi

while read sample; do
    cd ${sample}
        head -1 ${sample}.fastq | awk -v OFS='\t' -v s=$sample '{print substr($1, 2, length($1)-3), s}' >> ../sample_to_accesions_yeast
    cd ..
done < mapfile

output_file=2015_moller_scerevisiae_delcollection_notreatment_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_d

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}

cp ${output_file} ${output_file}.old

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' sample_to_accesions_yeast ${output_file}.old > ${output_file}

output_file=2015_moller_scerevisiae_delcollection_zeocin_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_z

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}

cp ${output_file} ${output_file}.old

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' sample_to_accesions_yeast ${output_file}.old > ${output_file}

output_file=2015_moller_scerevisiae_gap1circle_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_g

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}

cp ${output_file} ${output_file}.old

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' sample_to_accesions_yeast ${output_file}.old > ${output_file}

output_file=2015_moller_scerevisiae_wt_eccdnas.txt

if [ -f "${output_file}" ]; then
    rm ${output_file}
fi

while read sample; do
    cd ${sample}
        awk -v OFS='\t' -v s=${sample} '{print s, $1, $2, $3, $4}' ${sample}.ecc_caller_out.details.nolowq.txt >> ../${output_file}
    cd ..
done < mapfile_s

cp ${output_file} ${output_file}.old

sort -k2,2 -k3,3n ${output_file}.old > ${output_file}

cp ${output_file} ${output_file}.old

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' sample_to_accesions_yeast ${output_file}.old > ${output_file}