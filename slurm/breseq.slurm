#!/bin/bash
#SBATCH --job-name=breseq
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --account=fc_koskella
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=rrdebray@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/rrdebray/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/rrdebray/slurm_stderr/slurm-%j.out


WORK_DIR=
REFERENCE_GENOME=
THREADS=$SLURM_NTASKS
OUTPUT_FILE=
TRIMMED_READS1=
TRIMMED_READ2=

cd $WORK_DIR

$BRESEQ_PATH/breseq -r $REFERENCE_GENOME -j $THREADS -p -o $OUTPUT_FILE $TRIMMED_READS_1 $TRIMMED_READS2