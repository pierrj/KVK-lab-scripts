#!/bin/bash
#SBATCH --job-name=confirmsrs_sql_test_singledb
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/dummy_parallel

# ipcluster start -n 20 --cluster-id="cluster-id-dummy" &
# sleep 60
# ipython /global/home/users/pierrj/git/python/dummy_ipython_parallel.py
# ipcluster stop --cluster-id="cluster-id-dummy"

ipcontroller --nodb --cluster-id='cluster-id-dummy' &
sleep 60
for (( i = 0 ; i < 20; i++)); do ipengine --profile=pierrj --cluster-id='cluster-id-dummy' & done
sleep 60
ipython /global/home/users/pierrj/git/python/dummy_ipython_parallel.py
ipcluster stop --cluster-id="cluster-id-dummy"