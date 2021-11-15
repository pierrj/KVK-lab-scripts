#!/bin/bash
#SBATCH --job-name=ecc_sv_finder_moryzae
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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/structural_variation/moryzae

# list of moryzae genomes
MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/mapfile
# repbase + repeat modeler library on guy11 genome
RM_LIB=/global/scratch/users/pierrj/references/te_annotations/moryzae/repbase_plus_guy11_modeler_rm_database.fasta
# genome fasta files
GENOMES_PATH=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/
REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/guy11_genome_baoetal2017.fasta
# for mummer
PERCENT_ZEROES_FILTER=1

module load gnuplot
module load imagemagick
module load repeatmasker/4.1.1

# repeat mask all genomes, only for TEs not repeats
while read genome; do
    RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${GENOMES_PATH}/${genome}_genomic.fna
done < ${MAPFILE}

RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

if [ -f "jobqueue" ]; then
    rm jobqueue
fi

## run script that looks for svs that may have been caused by eccs, comparing guy11 vs all other genomes
while read genome;
do
    echo "/global/home/users/pierrj/git/bash/ecc_sv_finder.sh -r guy11_genome_baoetal2017.fasta.masked -q ${genome}_genomic.fna.masked -t 1 -n moryzae_v_${genome} -p ${PERCENT_ZEROES_FILTER} -o moryzae_v_${genome}_jpg -m moryzae_v_${genome}_tmp" >> jobqueue
    echo "/global/home/users/pierrj/git/bash/ecc_sv_finder.sh -r ${genome}_genomic.fna.masked -q guy11_genome_baoetal2017.fasta.masked -t 1 -n ${genome}_v_moryzae -p ${PERCENT_ZEROES_FILTER} -o ${genome}_v_moryzae_jpg -m ${genome}_v_moryzae_tmp" >> jobqueue
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