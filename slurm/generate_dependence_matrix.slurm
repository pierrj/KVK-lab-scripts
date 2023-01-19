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
#MIT License
#
#Copyright (c) 2023 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest

INPUT_DF=gene_info.full_model.rice_blast.txt

## rf params
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

## start by outputting set train and test data, in accordance with majority_fraction param
/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/output_train_test.py $INPUT_DF $MAJORITY_FRACTION $OUTPUT_COL_LIST $OUTPUT_XTRAIN $OUTPUT_XTEST

# fix line endings
tr -d '\015' <${OUTPUT_COL_LIST} > ${OUTPUT_COL_LIST}.fixed

if [ -f jobqueue_dep ]; then
    rm jobqueue_dep
fi

# parallel stuff
while read param; do
    echo "/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/dependencies_per_column.py $OUTPUT_XTRAIN $OUTPUT_XTEST $param $MAJORITY_FRACTION $APPROACH $ESTIMATORS $SPLIT $LEAF $FEATURES $DEPTH $BOOTSTRAP" >> jobqueue_dep
done < ${OUTPUT_COL_LIST}.fixed

if [ -f $OUTPUT_FILE ]; then
    rm $OUTPUT_FILE
fi

parallel -j ${SLURM_NTASKS} < jobqueue_dep >> ${OUTPUT_FILE}

## moving and reformating the files for plotting dependency matrix in R
mv $OUTPUT_FILE $OUTPUT_FILE.old

mv ${OUTPUT_COL_LIST}.fixed $OUTPUT_COL_LIST.fixed.old

cat <(echo x) <(echo f1) $OUTPUT_COL_LIST.fixed.old > $OUTPUT_COL_LIST

mv $OUTPUT_COL_LIST $OUTPUT_COL_LIST.old

tr "\n" "\t" < $OUTPUT_COL_LIST.old | sed 's/.$//' > header

cat header <(echo) $OUTPUT_FILE.old > $OUTPUT_FILE