#!/bin/bash
#SBATCH --job-name=eccdnavsrepeatcomparison_test
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/repeat_comparison_wrapper/

GENOME_FASTA=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
SAMPLE=G3
THREADS=${SLURM_NTASKS}
REPEAT_FILE=/global/scratch/users/pierrj/references/GUY11_pacbio_ET_KROJ.gff
ECCDNA_MAPFILE=G3.eccmapfile
SAMPLE_MAPFILE=G3.sample_mapfile

/global/home/users/pierrj/git/bash/ecc_dna_vs_repetitive.sh -g ${GENOME_FASTA} \
    -s ${SAMPLE} -t ${THREADS} -f ${REPEAT_FILE} \
    -e ${ECCDNA_MAPFILE} -m ${SAMPLE_MAPFILE}