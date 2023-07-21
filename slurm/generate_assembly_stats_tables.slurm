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


## wheat

GENOMES_PATH=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_to_annotate
MAPFILE=/global/scratch/users/pierrj/fungap_runs/wheat_blast/anne_curated_mapfile_annotated_busco.txt
host=wheat
acronym=MoT

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

busco -i genomes_for_busco -l sordariomycetes_odb10 -o ${host}.genomes_for_busco_buscofied -m genome -c $SLURM_NTASKS -f

# remove header
tail -n +2 ${host}.genomes_for_busco_buscofied/batch_summary.txt | sort -k1,1 > ${host}.busco.out

# process accession names just for wheat blast
awk -F "_" '{print $1$2}' $MAPFILE > accession_column

## assemble columns with paste and select useful columns
# from assembly stats
# N seqs
# N nucleotides
# N50
# L50
# from busco
# completeness
paste accession_column $(awk -F'\t' 'BEGIN{OFS="\t"} {print $3, $6, $21, $22}' genome_assembly_stats) $(awk -F'\t' 'BEGIN{OFS="\t"} {print $3}' ${host}.busco.out) > final_table

## add repeated values
awk -v host=$host -v acronym=$acronym -F'\t' 'BEGIN{OFS="\t"} {print $1, host, acronym, $2, $3, $4, $5, $6}' final_table > ${host}.final_table

cat <(echo -e isolate'\t'host'\t'acronym'\t'n_seqs'\t'n_bases'\t'N50'\t'L50'\t'busco_complete) \
    rice.final_table \
    wheat.final_table > pav_paper.assembly_stats.txt