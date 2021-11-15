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
#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.


cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

# permutation analysis of junction split reads for each eccdna category

## LTR ECCS

## grab all ltr eccdnas
group=ltr
group_file=ltr_eccs
file=ltr_eccdnas_all

if [ -f "${file}" ]; then
    rm ${file}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${group_file} >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas/${file}; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile


cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

## get eccdna counts (not junction split reads) and steps for rarefaction analysis
splitreadcount=$(wc -l ${file} | awk '{print $1}')
splitreadcount_steps=$(wc -l ${file} | awk '{print int($1/10)}')


for k in $(seq 0 99); do # 100 permutations
    if [ -f "${group}_regions_found_rarefaction_${k}" ]; then
        rm ${group}_regions_found_rarefaction_${k}
    fi
    for i in `seq ${splitreadcount_steps} ${splitreadcount_steps} ${splitreadcount}`; do ## 10 steps each
        shuf ${file} | head -n ${i} | sort -k1,1 -k2,2n | uniq | wc -l >> ${group}_regions_found_rarefaction_${k} ## sample according to steps
    done
done

# make big table
paste $(find . -maxdepth 1 -name "${group}_regions_found_rarefaction_*" | xargs -r ls -1 | tr "\n" " ") > ${group}_regions_found_rarefaction


## LARGE ECCS

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

## grab all large eccdnas
group=large
group_file=large_eccdnas
file=large_eccdnas_all

if [ -f "${file}" ]; then
    rm ${file}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${group_file} >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas/${file}; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile


cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

## get eccdna counts (not junction split reads) and steps for rarefaction analysis
splitreadcount=$(wc -l ${file} | awk '{print $1}')
splitreadcount_steps=$(wc -l ${file} | awk '{print int($1/10)}')



for k in $(seq 0 99); do # 100 permutations
    if [ -f "${group}_regions_found_rarefaction_${k}" ]; then
        rm ${group}_regions_found_rarefaction_${k}
    fi
    for i in `seq ${splitreadcount_steps} ${splitreadcount_steps} ${splitreadcount}`; do ## 10 steps each
        shuf ${file} | head -n ${i} | sort -k1,1 -k2,2n | uniq | wc -l >> ${group}_regions_found_rarefaction_${k} # sample according to steps
    done
done

# make big table
paste $(find . -maxdepth 1 -name "${group}_regions_found_rarefaction_*" | xargs -r ls -1 | tr "\n" " ") > ${group}_regions_found_rarefaction

## MICRO DNAS

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

## grab all micro dnas
group=micro
group_file=micro_dnas
file=micro_dnas_all

if [ -f "${file}" ]; then
    rm ${file}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${group_file} >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas/${file}; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile


cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/circrnas

## get eccdna counts (not junction split reads) and steps for rarefaction analysis
splitreadcount=$(wc -l ${file} | awk '{print $1}')
splitreadcount_steps=$(wc -l ${file} | awk '{print int($1/10)}')



for k in $(seq 0 99); do
    if [ -f "${group}_regions_found_rarefaction_${k}" ]; then
        rm ${group}_regions_found_rarefaction_${k}
    fi
    for i in `seq ${splitreadcount_steps} ${splitreadcount_steps} ${splitreadcount}`; do ## 10 steps each
        shuf ${file} | head -n ${i} | sort -k1,1 -k2,2n | uniq | wc -l >> ${group}_regions_found_rarefaction_${k} # sample according to steps
    done
done

paste $(find . -maxdepth 1 -name "${group}_regions_found_rarefaction_*" | xargs -r ls -1 | tr "\n" " ") > ${group}_regions_found_rarefaction