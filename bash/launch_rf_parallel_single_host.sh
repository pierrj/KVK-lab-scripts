## rice full model ##

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest

INPUT_DF=gene_info.full_model.rice_blast.txt
OUTPUT_FILE=rf_results_replicated.${INPUT_DF}

MAJORITY_FRACTION=0.5
APPROACH=RF
ESTIMATORS=2000
SPLIT=2
LEAF=1
FEATURES=None
DEPTH=None
BOOTSTRAP=True

if [ -f $OUTPUT_FILE ]; then
    rm $OUTPUT_FILE
fi

if [ -f jobqueue ]; then
    rm jobqueue
fi

REPLICATES=100

for i in $(seq $REPLICATES); do
    echo "/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_perf_test_parallel.py $INPUT_DF $MAJORITY_FRACTION $APPROACH $ESTIMATORS $SPLIT $LEAF $FEATURES $DEPTH $BOOTSTRAP" >> jobqueue
done

N_NODES=3

split --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

for node in $(seq -f "%02g" 1 ${N_NODES})
do
    sbatch --job-name=$node.rf --export=ALL,OUTPUT_FILE=$OUTPUT_FILE,node=$node \
    --account=co_minium \
    /global/home/users/pierrj/git/slurm/htc4_gnu_parallel_rf.slurm
done

mv $OUTPUT_FILE ${OUTPUT_FILE}.old

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/average_rf_results.py ${OUTPUT_FILE}.old $OUTPUT_FILE


## rice reduced model ##

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/cross_host

INPUT_DF=gene_info.cross_host.rice_blast.txt
OUTPUT_FILE=rf_results_replicated.${INPUT_DF}

MAJORITY_FRACTION=0.5
APPROACH=RF
ESTIMATORS=2000
SPLIT=2
LEAF=1
FEATURES=None
DEPTH=None
BOOTSTRAP=True

if [ -f $OUTPUT_FILE ]; then
    rm $OUTPUT_FILE
fi

if [ -f jobqueue ]; then
    rm jobqueue
fi

REPLICATES=100

for i in $(seq $REPLICATES); do
    echo "/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_perf_test_parallel.py $INPUT_DF $MAJORITY_FRACTION $APPROACH $ESTIMATORS $SPLIT $LEAF $FEATURES $DEPTH $BOOTSTRAP" >> jobqueue
done

N_NODES=3

split --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

for node in $(seq -f "%02g" 1 ${N_NODES})
do
    sbatch --job-name=$node.rf --export=ALL,OUTPUT_FILE=$OUTPUT_FILE,node=$node \
    --account=co_minium \
    /global/home/users/pierrj/git/slurm/htc4_gnu_parallel_rf.slurm
done

mv $OUTPUT_FILE ${OUTPUT_FILE}.old

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/average_rf_results.py ${OUTPUT_FILE}.old $OUTPUT_FILE


## wheat reduced model ##

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/cross_host

INPUT_DF=gene_info.cross_host.wheat_blast.txt
OUTPUT_FILE=rf_results_replicated.${INPUT_DF}

MAJORITY_FRACTION=0.5
APPROACH=RF
ESTIMATORS=2000
SPLIT=2
LEAF=1
FEATURES=None
DEPTH=None
BOOTSTRAP=True

if [ -f $OUTPUT_FILE ]; then
    rm $OUTPUT_FILE
fi

if [ -f jobqueue ]; then
    rm jobqueue
fi

REPLICATES=100

for i in $(seq $REPLICATES); do
    echo "/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_perf_test_parallel.py $INPUT_DF $MAJORITY_FRACTION $APPROACH $ESTIMATORS $SPLIT $LEAF $FEATURES $DEPTH $BOOTSTRAP" >> jobqueue
done

N_NODES=3

split --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

for node in $(seq -f "%02g" 1 ${N_NODES})
do
    sbatch --job-name=$node.rf --export=ALL,OUTPUT_FILE=$OUTPUT_FILE,node=$node \
    --account=co_minium \
    /global/home/users/pierrj/git/slurm/htc4_gnu_parallel_rf.slurm
done

mv $OUTPUT_FILE ${OUTPUT_FILE}.old

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/average_rf_results.py ${OUTPUT_FILE}.old $OUTPUT_FILE
