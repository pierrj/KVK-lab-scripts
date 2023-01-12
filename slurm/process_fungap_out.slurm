#!/bin/bash
#SBATCH --job-name=process_fungap_out
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
#MIT License
#
#Copyright (c) 2023 Pierre Michel Joubert
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

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap

if [ ! -d "all_proteomes" ]; then
    mkdir all_proteomes
fi

MAPFILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/genomes_mapfile
FUNGAP_DIR=/global/scratch/users/pierrj/fungap_runs/gladieux_all/

# copy fungap out files over
while read genome; do
    cp ${FUNGAP_DIR}/${genome}/fungap_out/fungap_out/fungap_out_prot.faa all_proteomes/${genome}_fungap_out_prot.faa
done < $MAPFILE

if [ ! -d "all_proteomes_corrected" ]; then
    mkdir all_proteomes_corrected
fi

source activate /global/scratch/users/pierrj/conda_envs/orthofinder
module purge

## process proteome sequences for use w/ orthofinder
/global/scratch/users/pierrj/conda_envs/orthofinder/bin/python /global/home/users/pierrj/git/python/process_protein_sequences_for_orthofinder.py all_proteomes \
    all_proteomes_corrected \
    gladieux_et_al_assembly_stats.txt \
    output_accession_names

cd /global/scratch/users/pierrj/fungap_runs/wheat_blast

MAPFILE=/global/scratch/users/pierrj/fungap_runs/wheat_blast/wheat_blast_busco_greater_than_90_annotated

## some issues with the gff file that need to be fixed before processing
while read genome; do
    grep -v '#' ${genome}/fungap_out/fungap_out/fungap_out.gff3 > ${genome}/fungap_out/fungap_out/fungap_out.nocomments.gff3
    sed -i '1s/^/##gff-version 3\n/' ${genome}/fungap_out/fungap_out/fungap_out.nocomments.gff3
done < $MAPFILE

# process gffs for use w/ orthofinder
/global/scratch/users/pierrj/conda_envs/orthofinder/bin/python /global/home/users/pierrj/git/python/process_gffs_for_orthofinder.py $MAPFILE

# fix gffs file so they are gff3 formatted
while read genome; do
    mv all_gffs_fixed/${genome}_fungap_out.fixed.gff3 all_gffs_fixed/${genome}_fungap_out.fixed.gff3.old
    grep -v gff-version=3=true all_gffs_fixed/${genome}_fungap_out.fixed.gff3.old > all_gffs_fixed/${genome}_fungap_out.fixed.gff3
    rm all_gffs_fixed/${genome}_fungap_out.fixed.gff3.old
done < $MAPFILE