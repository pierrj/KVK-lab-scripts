
cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/genome_exclusion_test

INPUT_DF=../gene_info.genes_tes_gc.lengths_effectors.guy11_histone_expression.go_pfam_eccdnas.txt
OUTPUT_FILE=rf_results.genome_exclusion_test
MAJORITY_FRACTION=1.0
APPROACH=SMOTE
ESTIMATORS=2000
SPLIT=2
LEAF=1
FEATURES=sqrt
DEPTH=None
BOOTSTRAP=True

if [ -f $OUTPUT_FILE ]; then
    rm $OUTPUT_FILE
fi

if [ -f jobqueue ]; then
    rm jobqueue
fi

REPLICATES=10

for i in $(seq $REPLICATES); do
    echo "/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_perf_test_genome_exclusion.py $INPUT_DF $MAJORITY_FRACTION $APPROACH $ESTIMATORS $SPLIT $LEAF $FEATURES $DEPTH $BOOTSTRAP $INPUT_DF2" >> jobqueue
done

CPUS=10

split --number=l/${CPUS} --numeric-suffixes=1 jobqueue jobqueue_

for cpu in $(seq -f "%02g" 1 ${CPUS})
do
    # echo $cpu
    sbatch --job-name=$cpu.rf -n 1 --export=ALL,OUTPUT_FILE=$OUTPUT_FILE,cpu=$cpu \
    --account=fc_kvkallow -p savio2_htc --qos=savio_normal \
    /global/home/users/pierrj/git/slurm/htc4_gnu_parallel_rf_single_cpu.slurm
done

## after everything is done
mv $OUTPUT_FILE ${OUTPUT_FILE}.old

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/average_rf_results.py ${OUTPUT_FILE}.old $OUTPUT_FILE






## GRID SEARCH ##

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest

# print param grid #
INPUT_DF=gene_info.genes_tes_gc.lengths_effectors.guy11_histone_expression.go_pfam_eccdnas.txt
OUTPUT_FILE=test7.rf_results_grid_search.${INPUT_DF}

if [ -f $OUTPUT_FILE ]; then
    rm ${OUTPUT_FILE}*
fi

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/print_param_grid.py > param_grid

shuf param_grid > param_grid_shuffled

cat param_grid_shuffled param_grid_shuffled param_grid_shuffled param_grid_shuffled param_grid_shuffled param_grid_shuffled > param_grid_shuffled_replicates

if [ -f jobqueue ]; then
    rm jobqueue*
fi

while read APPROACH MAJORITY_FRACTION ESTIMATORS SPLIT LEAF FEATURES DEPTH BOOTSTRAP; do
    echo "/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_perf_test_parallel.py $INPUT_DF $MAJORITY_FRACTION $APPROACH $ESTIMATORS $SPLIT $LEAF $FEATURES $DEPTH $BOOTSTRAP" >> jobqueue
done < param_grid_shuffled_replicates

N_NODES=4

split --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

for node in $(seq -f "%02g" 1 ${N_NODES})
do
    sbatch --job-name=$node.rf --export=ALL,OUTPUT_FILE=$OUTPUT_FILE,node=$node \
    --account=co_minium \
    /global/home/users/pierrj/git/slurm/htc4_gnu_parallel_rf.slurm
done


CPUS=181

split --number=l/${CPUS} --numeric-suffixes=1 jobqueue jobqueue_

for cpu in $(seq -f "%03g" 1 ${CPUS})
do
    # echo $cpu
    sbatch --job-name=$cpu.rf --mem 8000M -n 1 --export=ALL,OUTPUT_FILE=$OUTPUT_FILE,cpu=$cpu \
    --account=co_minium \
    /global/home/users/pierrj/git/slurm/htc4_gnu_parallel_rf_single_cpu.slurm
done

N_NODES=31

split --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

for node in $(seq -f "%02g" 1 ${N_NODES})
do
    sbatch --job-name=$node.rf --requeue --export=ALL,OUTPUT_FILE=$OUTPUT_FILE,node=$node \
    --account=co_minium \
    /global/home/users/pierrj/git/slurm/gnu_parallel_rf_low_prio.slurm
done



for node in $(seq -f "%02g" 1 ${N_NODES})
do
    wc -l jobqueue_${node}
done

## after everything is done

cat ${OUTPUT_FILE}_* > ${OUTPUT_FILE}

mv $OUTPUT_FILE ${OUTPUT_FILE}.old

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/average_rf_results.py ${OUTPUT_FILE}.old $OUTPUT_FILE
