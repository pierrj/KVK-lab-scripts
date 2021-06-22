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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/structural_variation/moryzae

MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/mapfile
RM_LIB=/global/scratch/users/pierrj/references/te_annotations/moryzae/repbase_plus_guy11_modeler_rm_database.fasta
GENOMES_PATH=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/
REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/guy11_genome_baoetal2017.fasta
PERCENT_ZEROES_FILTER=1

module load gnuplot
module load imagemagick
module load repeatmasker/4.1.1

while read genome; do
    RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${GENOMES_PATH}/${genome}_genomic.fna
done < ${MAPFILE}

RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

if [ -f "jobqueue" ]; then
    rm jobqueue
fi

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

while read output; do
cd ${output}
    cp *.pdf ../output_pdfs
cd ..
done < output_mapfile