
cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest

INPUT_DF=gene_info.genes_tes_gc.lengths_effectors.guy11_histone_expression.txt
OUTPUT_FILE=rf_results.genes_tes_gc.lengths_effectors.guy11_histone_expression.txt

if [ -f $OUTPUT_FILE ]; then
    rm $OUTPUT_FILE
fi

if [ -f jobqueue ]; then
    rm jobqueue
fi

REPLICATES=10
for APPROACH in "RF" "SMOTE" "BRFC" "RF_balanced" "RF_balanced_subsample"; do
    for MAJORITY_FRACTION in 0.05 0.1 0.25 0.5 0.75 0.95; do
        for i in $(seq $REPLICATES); do
            echo "/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_perf_test_parallel.py $INPUT_DF $MAJORITY_FRACTION $APPROACH" >> jobqueue
        done
    done
done

mv jobqueue jobqueue_old

shuf jobqueue_old > jobqueue

cores=56
mem=224000M
N_NODES=4

split --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

for node in $(seq -f "%02g" 1 ${N_NODES})
do
    sbatch --job-name=$node.rf --mem=${mem} -n $cores --export=ALL,OUTPUT_FILE=$OUTPUT_FILE,node=$node \
    -p savio4_htc --qos=minium_htc4_normal --account=co_minium \
    /global/home/users/pierrj/git/slurm/htc4_gnu_parallel_rf.slurm
done

## after everything is done
mv $OUTPUT_FILE ${OUTPUT_FILE}.old

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/average_rf_results.py ${OUTPUT_FILE}.old $OUTPUT_FILE