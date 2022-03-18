#!/bin/bash
#SBATCH --job-name=ecc_sv_finder_mphaseolina
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/scooping_emile/

# list of moryzae genomes
MAPFILE=/global/scratch/users/pierrj/scooping_emile/genomes_mapfile
# repbase + repeat modeler library on guy11 genome
RM_LIB=/global/scratch/users/pierrj/repbase/RepBase25.04.fasta/fngrep.ref
# genome fasta files
GENOMES_PATH=/global/scratch/users/pierrj/scooping_emile/genomes
# for mummer
PERCENT_ZEROES_FILTER=1

# # for download
# GENOME_LINKS=/global/scratch/users/pierrj/scooping_emile/genomes_links
# cd ${GENOMES_PATH}

# while read genome; do
#     wget ${genome}
# done < ${GENOME_LINKS}

# gunzip *

# cd /global/scratch/users/pierrj/scooping_emile/

module load gnuplot
module load imagemagick
module load repeatmasker/4.1.1


# repeat mask all genomes, only for TEs not repeats
while read genome; do
    RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${GENOMES_PATH}/${genome}_genomic.fna
done < ${MAPFILE}

if [ -f "jobqueue" ]; then
    rm jobqueue
fi

## run script that looks for svs that may have been caused by eccs, comparing guy11 vs all other genomes
## do all by all
while read reference;
do
    while read genome;
    do
        echo "/global/home/users/pierrj/git/bash/ecc_sv_finder.sh -r ${reference}_genomic.fna.masked -q ${genome}_genomic.fna.masked -t 1 -n ${reference}_v_${genome} -p ${PERCENT_ZEROES_FILTER} -o ${reference}_v_${genome}_jpg -m ${reference}_v_${genome}_tmp" >> jobqueue
        echo "/global/home/users/pierrj/git/bash/ecc_sv_finder.sh -r ${genome}_genomic.fna.masked -q ${reference}_genomic.fna.masked -t 1 -n ${genome}_v_${reference} -p ${PERCENT_ZEROES_FILTER} -o ${genome}_v_${reference}_jpg -m ${genome}_v_${reference}_tmp" >> jobqueue
    done < ${MAPFILE}
done < ${MAPFILE}

parallel -j ${SLURM_NTASKS} < jobqueue

if [ -d "output_pdfs" ]; then
    rm -r output_pdfs
fi

mkdir output_pdfs

ls | grep jpg > output_mapfile

## move files around
while read output; do
cd ${output}
    cp *.pdf ../output_pdfs
cd ..
done < output_mapfile