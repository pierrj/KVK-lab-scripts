#!/bin/bash
#SBATCH --job-name=permutation_eccdna_forming_regions_by_type
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

## LTR ECCS

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

group=ltr
group_file=ltr_eccs
file=ltr_eccdnas_all

if [ -f "${file}" ]; then
    rm ${file}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${group_file} >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas/${file}; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile


cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

splitreadcount=$(wc -l ${file} | awk '{print $1}')
splitreadcount_steps=$(wc -l ${file} | awk '{print int($1/10)}')



for k in $(seq 0 99); do
    if [ -f "${group}_regions_found_rarefaction_${k}" ]; then
        rm ${group}_regions_found_rarefaction_${k}
    fi
    for i in `seq ${splitreadcount_steps} ${splitreadcount_steps} ${splitreadcount}`; do
        shuf ${file} | head -n ${i} | sort -k1,1 -k2,2n | uniq | wc -l >> ${group}_regions_found_rarefaction_${k}
    done
done

paste $(find . -maxdepth 1 -name "${group}_regions_found_rarefaction_*" | xargs -r ls -1 | tr "\n" " ") > ${group}_regions_found_rarefaction


## LARGE ECCS

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

group=large
group_file=large_eccdnas
file=large_eccdnas_all

if [ -f "${file}" ]; then
    rm ${file}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${group_file} >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas/${file}; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile


cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

splitreadcount=$(wc -l ${file} | awk '{print $1}')
splitreadcount_steps=$(wc -l ${file} | awk '{print int($1/10)}')



for k in $(seq 0 99); do
    if [ -f "${group}_regions_found_rarefaction_${k}" ]; then
        rm ${group}_regions_found_rarefaction_${k}
    fi
    for i in `seq ${splitreadcount_steps} ${splitreadcount_steps} ${splitreadcount}`; do
        shuf ${file} | head -n ${i} | sort -k1,1 -k2,2n | uniq | wc -l >> ${group}_regions_found_rarefaction_${k}
    done
done

paste $(find . -maxdepth 1 -name "${group}_regions_found_rarefaction_*" | xargs -r ls -1 | tr "\n" " ") > ${group}_regions_found_rarefaction

## MICRO DNAS

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

group=micro
group_file=micro_dnas
file=micro_dnas_all

if [ -f "${file}" ]; then
    rm ${file}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${group_file} >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas/${file}; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile


cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

splitreadcount=$(wc -l ${file} | awk '{print $1}')
splitreadcount_steps=$(wc -l ${file} | awk '{print int($1/10)}')



for k in $(seq 0 99); do
    if [ -f "${group}_regions_found_rarefaction_${k}" ]; then
        rm ${group}_regions_found_rarefaction_${k}
    fi
    for i in `seq ${splitreadcount_steps} ${splitreadcount_steps} ${splitreadcount}`; do
        shuf ${file} | head -n ${i} | sort -k1,1 -k2,2n | uniq | wc -l >> ${group}_regions_found_rarefaction_${k}
    done
done

paste $(find . -maxdepth 1 -name "${group}_regions_found_rarefaction_*" | xargs -r ls -1 | tr "\n" " ") > ${group}_regions_found_rarefaction