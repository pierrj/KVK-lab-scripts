#!/bin/bash
#SBATCH --job-name=repeatmask_per_genome
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

DATA_OUTPUT_DIR=repeatmasker

cd $PROJ_DIR

if [ -d $DATA_OUTPUT_DIR ]; then
    rm -r $DATA_OUTPUT_DIR
fi

mkdir $DATA_OUTPUT_DIR

cd $DATA_OUTPUT_DIR

module load repeatmasker/4.1.1

while read genome; do
    RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${GENOME_FASTA}
done < ${MAPFILE}

while read genome; do
    awk '{print $11}' ${GENOME_FASTA}.out | paste ${GENOME_FASTA}.out.gff - > ${GENOME_FASTA}.w_class.gff
done < ${MAPFILE}