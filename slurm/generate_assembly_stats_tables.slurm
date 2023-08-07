#!/bin/bash
#SBATCH --job-name=get_assembly_stats
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

source activate /global/scratch/users/pierrj/conda_envs/gaas

export PERL5LIB=/global/scratch/users/pierrj/conda_envs/gaas/lib

cd /global/scratch/users/pierrj/PAV_SV/PAV/assembly_stats

## rice first

GENOMES_PATH=/global/scratch/users/pierrj/fungap_runs/gladieux_all/genomes_to_annotate/
MAPFILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/genomes_mapfile_no_mgrisae
host=rice
acronym=MoO

# assembly stats
while read genome; do
    gaas_fasta_statistics.pl -f ${GENOMES_PATH}/${genome}_genomic.fna | awk '{print $(NF-1)}' | awk '$1 !~ /--/' > ${genome}.stats
done < ${MAPFILE}

if [ -f "genome_assembly_stats" ]; then
    rm genome_assembly_stats
fi

while read genome; do
    tr '\n' '\t' < ${genome}.stats >> genome_assembly_stats
    echo >> genome_assembly_stats
done < ${MAPFILE}

# busco
if [ -d "genomes_for_busco" ]; then
    rm -r genomes_for_busco
fi

mkdir genomes_for_busco

while read genome; do
    cp ${GENOMES_PATH}/${genome}_genomic.fna genomes_for_busco/
done < ${MAPFILE}

source activate /global/scratch/users/pierrj/conda_envs/BUSCO

export PATH=/global/scratch/users/pierrj/conda_envs/BUSCO/bin:$PATH

# busco -i genomes_for_busco -l sordariomycetes_odb10 -o ${host}.genomes_for_busco_buscofied -m genome -c $SLURM_NTASKS -f

# remove header
tail -n +2 ${host}.genomes_for_busco_buscofied/batch_summary.txt | sort -k1,1 > ${host}.busco.out

## assemble columns with paste and select useful columns
# from assembly stats
# N seqs
# N nucleotides
# N50
# L50
# from busco
# completeness
paste $MAPFILE <(awk -F'\t' 'BEGIN{OFS="\t"} {print $3, $6, $21, $22}' genome_assembly_stats) <(awk -F'\t' 'BEGIN{OFS="\t"} {print $3}' ${host}.busco.out) > final_table

## add repeated values
awk -v host=$host -v acronym=$acronym -F'\t' 'BEGIN{OFS="\t"} {print $1, host, acronym, $2, $3, $4, $5, $6}' final_table > ${host}.final_table


source activate /global/scratch/users/pierrj/conda_envs/gaas

export PERL5LIB=/global/scratch/users/pierrj/conda_envs/gaas/lib

cd /global/scratch/users/pierrj/PAV_SV/PAV/assembly_stats

## wheat

GENOMES_PATH=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_to_annotate
MAPFILE=/global/scratch/users/pierrj/fungap_runs/wheat_blast/anne_curated_mapfile_annotated_busco.txt
host=wheat
acronym=MoT

# sort mapfile
sort -k1,1 $MAPFILE > wheat_mapfile
MAPFILE=wheat_mapfile

# assembly stats
while read genome; do
    gaas_fasta_statistics.pl -f ${GENOMES_PATH}/${genome}_genomic.fna | awk '{print $(NF-1)}' | awk '$1 !~ /--/' > ${genome}.stats
done < ${MAPFILE}

if [ -f "genome_assembly_stats" ]; then
    rm genome_assembly_stats
fi

while read genome; do
    tr '\n' '\t' < ${genome}.stats >> genome_assembly_stats
    echo >> genome_assembly_stats
done < ${MAPFILE}

# busco
if [ -d "genomes_for_busco" ]; then
    rm -r genomes_for_busco
fi

mkdir genomes_for_busco

while read genome; do
    cp ${GENOMES_PATH}/${genome}_genomic.fna genomes_for_busco/
done < ${MAPFILE}

source activate /global/scratch/users/pierrj/conda_envs/BUSCO

export PATH=/global/scratch/users/pierrj/conda_envs/BUSCO/bin:$PATH

# busco -i genomes_for_busco -l sordariomycetes_odb10 -o ${host}.genomes_for_busco_buscofied -m genome -c $SLURM_NTASKS -f

# remove header
tail -n +2 ${host}.genomes_for_busco_buscofied/batch_summary.txt | sort -k1,1 > ${host}.busco.out

## assemble columns with paste and select useful columns
# from assembly stats
# N seqs
# N nucleotides
# N50
# L50
# from busco
# completeness

awk -F "_" '{print $1$2}' $MAPFILE > wheat_colnames
MAPFILE=wheat_colnames

paste $MAPFILE <(awk -F'\t' 'BEGIN{OFS="\t"} {print $3, $6, $21, $22}' genome_assembly_stats) <(awk -F'\t' 'BEGIN{OFS="\t"} {print $3}' ${host}.busco.out) > final_table

## add repeated values
awk -v host=$host -v acronym=$acronym -F'\t' 'BEGIN{OFS="\t"} {print $1, host, acronym, $2, $3, $4, $5, $6}' final_table > ${host}.final_table

cat <(echo -e isolate'\t'host'\t'acronym'\t'n_seqs'\t'n_bases'\t'N50'\t'L50'\t'busco_complete) \
    rice.final_table \
    wheat.final_table > pav_paper.assembly_stats.txt