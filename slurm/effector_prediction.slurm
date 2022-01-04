#!/bin/bash
#SBATCH --job-name=effector_prediction
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/effector_prediction


INPUT_FILE=guy11_fungap_out_12_28_20_prot.faa
OUTPUT_NAME=guy11

module load seqtk

cat ${INPUT_FILE} | /global/scratch/users/pierrj/tmhmm-2.0c/bin/tmhmm > ${OUTPUT_NAME}.tmhmmout

awk '{if ($7 == 0){print $2}} ' ${OUTPUT_NAME}.tmhmmout > ${OUTPUT_NAME}.notm.names

seqtk subseq ${INPUT_FILE} ${OUTPUT_NAME}.notm.names > ${OUTPUT_NAME}.notm.faa

/global/scratch/users/pierrj/signalp-4.1/signalp -t euk -u 0.34 -U 0.34 -f short ${OUTPUT_NAME}.notm.faa > ${OUTPUT_NAME}.signalpout.short

awk '{ if ($10 == "Y") {print $1} }' ${OUTPUT_NAME}.signalpout.short > ${OUTPUT_NAME}.sp.names

seqtk subseq ${OUTPUT_NAME}.notm.faa ${OUTPUT_NAME}.sp.names > ${OUTPUT_NAME}.notm.sp.faa

python /global/scratch/users/pierrj/EffectorP-3.0/EffectorP.py -i ${OUTPUT_NAME}.notm.sp.faa > ${OUTPUT_NAME}.effectorpout

awk '{ if ($(NF) == "effector") {print $1}}' ${OUTPUT_NAME}.effectorpout > ${OUTPUT_NAME}.predicted_effectors