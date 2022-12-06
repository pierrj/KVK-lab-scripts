## FULL MODEL
cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest

INPUT_DF=gene_info.full_model.rice_blast.txt
OUTPUT_FILE=rf_importances_replicated.${INPUT_DF}

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
    echo "/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_importances_parallel.py $INPUT_DF $MAJORITY_FRACTION $APPROACH $ESTIMATORS $SPLIT $LEAF $FEATURES $DEPTH $BOOTSTRAP" >> jobqueue
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

grep -v "any_te" ${OUTPUT_FILE}.old > ${OUTPUT_FILE}.old.nocolnames

grep "any_te" ${OUTPUT_FILE}.old | sort | uniq > ${OUTPUT_FILE}.old.colnames

cat ${OUTPUT_FILE}.old.colnames ${OUTPUT_FILE}.old.nocolnames > ${OUTPUT_FILE}.old

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/average_importances_results.py ${OUTPUT_FILE}.old $OUTPUT_FILE

mv $OUTPUT_FILE ${OUTPUT_FILE}.old

## transpose
awk -v OFS=';' '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' $OUTPUT_FILE.old > $OUTPUT_FILE.transposed

cat $(echo -e Feature'\t'Importance ) $OUTPUT_FILE.transposed > $OUTPUT_FILE

mv $OUTPUT_FILE $OUTPUT_FILE.old

tr [:blank:] \\t <$OUTPUT_FILE.old > $OUTPUT_FILE



## REDUCED MODEL
cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/cross_host

INPUT_DF=gene_info.cross_host.rice_blast.txt
OUTPUT_FILE=rf_importances_replicated.${INPUT_DF}

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
    echo "/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_importances_parallel.py $INPUT_DF $MAJORITY_FRACTION $APPROACH $ESTIMATORS $SPLIT $LEAF $FEATURES $DEPTH $BOOTSTRAP" >> jobqueue
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

grep -v "any_te" ${OUTPUT_FILE}.old > ${OUTPUT_FILE}.old.nocolnames

grep "any_te" ${OUTPUT_FILE}.old | sort | uniq > ${OUTPUT_FILE}.old.colnames

cat ${OUTPUT_FILE}.old.colnames ${OUTPUT_FILE}.old.nocolnames > ${OUTPUT_FILE}.old

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/average_importances_results.py ${OUTPUT_FILE}.old $OUTPUT_FILE

mv $OUTPUT_FILE ${OUTPUT_FILE}.old

## transpose
awk -v OFS='\t' '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' $OUTPUT_FILE.old > $OUTPUT_FILE.transposed

cat <(echo -e Feature'\t'Importance ) $OUTPUT_FILE.transposed > $OUTPUT_FILE

mv $OUTPUT_FILE $OUTPUT_FILE.old

tr [:blank:] \\t <$OUTPUT_FILE.old > $OUTPUT_FILE



## WHEAT BLAST MODEL
cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/cross_host

INPUT_DF=gene_info.cross_host.wheat_blast.txt
OUTPUT_FILE=rf_importances_replicated.${INPUT_DF}

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
    echo "/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/rf_importances_parallel.py $INPUT_DF $MAJORITY_FRACTION $APPROACH $ESTIMATORS $SPLIT $LEAF $FEATURES $DEPTH $BOOTSTRAP" >> jobqueue
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

grep -v "any_te" ${OUTPUT_FILE}.old > ${OUTPUT_FILE}.old.nocolnames

grep "any_te" ${OUTPUT_FILE}.old | sort | uniq > ${OUTPUT_FILE}.old.colnames

cat ${OUTPUT_FILE}.old.colnames ${OUTPUT_FILE}.old.nocolnames > ${OUTPUT_FILE}.old

/global/scratch/users/pierrj/conda_envs/random_forest/bin/python /global/home/users/pierrj/git/python/average_importances_results.py ${OUTPUT_FILE}.old $OUTPUT_FILE

mv $OUTPUT_FILE ${OUTPUT_FILE}.old

## transpose
awk -v OFS='\t' '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' $OUTPUT_FILE.old > $OUTPUT_FILE.transposed

cat <(echo -e Feature'\t'Importance ) $OUTPUT_FILE.transposed > $OUTPUT_FILE

mv $OUTPUT_FILE $OUTPUT_FILE.old

tr [:blank:] \\t <$OUTPUT_FILE.old > $OUTPUT_FILE