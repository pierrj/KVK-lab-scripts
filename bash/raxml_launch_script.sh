#!/bin/bash


cd /global/scratch/users/pierrj/PAV_SV/PAV/raxml_ng_test

if [ -f "jobqueue_search" ]; then
    rm jobqueue_search
fi

echo -e search_3'\t'33333 >> jobqueue_search
echo -e search_4'\t'44444 >> jobqueue_search
echo -e search_5'\t'55555 >> jobqueue_search

while read -r jobname seed; do
    export jobname=$jobname
    export seed=$seed
    envsubst < /global/home/users/pierrj/git/slurm/raxml_search_parallel.slurm > ${jobname}.slurm
    sbatch ${jobname}.slurm
done < jobqueue_search

if [ -f "jobqueue_bs" ]; then
    rm jobqueue_bs
fi

echo -e bs_1'\t'11111 >> jobqueue_bs
echo -e bs_2'\t'22222 >> jobqueue_bs
echo -e bs_3'\t'33333 >> jobqueue_bs
echo -e bs_4'\t'44444 >> jobqueue_bs
echo -e bs_5'\t'55555 >> jobqueue_bs


while read -r jobname seed; do
    export jobname=$jobname
    export seed=$seed
    envsubst < /global/home/users/pierrj/git/slurm/raxml_bootstrapping_parallel.slurm > ${jobname}.slurm
    sbatch ${jobname}.slurm
done < jobqueue_bs

