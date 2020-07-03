#!/bin/bash
#SBATCH --job-name=test_memory
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/anygenome_rice_test
SAMPLE="RC_1A"
THREADS="24"

ipcontroller --profile=pierrj --nodb --cluster-id='cluster-id-memory-test' &
sleep 60
for (( i = 0 ; i < ${THREADS}; i++)); do ipengine --profile=pierrj --cluster-id='cluster-id-memory-test' & done
sleep 60
ipython /global/home/users/pierrj/git/python/file_opening_memory_test.py lengthfiltered.merged.splitreads.${SAMPLE}.renamed.bed outwardfacing.${SAMPLE}.renamed.bed ${SAMPLE} 57
ipcluster stop --cluster-id="cluster-id-memory-test" --profile=pierrj


# ipcluster start -n 4 --cluster-id="cluster-id-memory-test" --profile=pierrj &
# sleep 60
# ipython /global/home/users/pierrj/git/python/file_opening_memory_test.py lengthfiltered.merged.splitreads.${SAMPLE}.renamed.bed outwardfacing.${SAMPLE}.renamed.bed ${SAMPLE} 56
# ipcluster stop --cluster-id="cluster-id-memory-test" --profile=pierrj