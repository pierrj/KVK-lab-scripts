#!/bin/bash
#SBATCH --job-name=generate_dependence_matrix
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --account=co_minium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest

INPUT_DF=gene_info.full_model.rice_blast.txt

MAJORITY_FRACTION=0.5
APPROACH=RF
ESTIMATORS=2000
SPLIT=2
LEAF=1
FEATURES=None
DEPTH=None
BOOTSTRAP=True

OUTPUT_COL_LIST=parameters_list.txt
OUTPUT_XTRAIN=x_train.${INPUT_DF}
OUTPUT_XTEST=x_test.${INPUT_DF}

OUTPUT_FILE=dependency_matrix.${INPUT_DF}

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/output_train_test.py $INPUT_DF $MAJORITY_FRACTION $OUTPUT_COL_LIST $OUTPUT_XTRAIN $OUTPUT_XTEST


tr -d '\015' <${OUTPUT_COL_LIST} > ${OUTPUT_COL_LIST}.fixed

if [ -f jobqueue_dep ]; then
    rm jobqueue_dep
fi

while read param; do
    echo "/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/dependencies_per_column.py $OUTPUT_XTRAIN $OUTPUT_XTEST $param $MAJORITY_FRACTION $APPROACH $ESTIMATORS $SPLIT $LEAF $FEATURES $DEPTH $BOOTSTRAP" >> jobqueue_dep
done < ${OUTPUT_COL_LIST}.fixed

if [ -f $OUTPUT_FILE ]; then
    rm $OUTPUT_FILE
fi

parallel -j ${SLURM_NTASKS} < jobqueue_dep >> ${OUTPUT_FILE}

mv $OUTPUT_FILE $OUTPUT_FILE.old

mv ${OUTPUT_COL_LIST}.fixed $OUTPUT_COL_LIST.fixed.old

cat <(echo x) <(echo f1) $OUTPUT_COL_LIST.fixed.old > $OUTPUT_COL_LIST

mv $OUTPUT_COL_LIST $OUTPUT_COL_LIST.old

tr "\n" "\t" < $OUTPUT_COL_LIST.old | sed 's/.$//' > header

cat header <(echo) $OUTPUT_FILE.old > $OUTPUT_FILE