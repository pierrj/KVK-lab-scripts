

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/go

PROTEOMES_PATH=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/all_proteomes_corrected/
MAPFILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/proteomes_mapfile_no_mgrisae


if [ -f jobqueue ]; then
    rm jobqueue
fi

while read proteome; do
    echo "/global/scratch/users/pierrj/conda_envs/pannzer/bin/python /global/scratch/users/pierrj/pannzer2/SANSPANZ.3/runsanspanz.py -R -o ',DE.${proteome}.out,GO.${proteome}.out,anno.${proteome}.out' -s 'Pyricularia oryzae' < ${PROTEOMES_PATH}/${proteome}" >> jobqueue
done < $MAPFILE

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
    /global/home/users/pierrj/git/slurm/htc4_gnu_parallel_pannzer.slurm
done


cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/go

## actually just do single cores

PROTEOMES_PATH=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/all_proteomes_corrected/
MAPFILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/proteomes_mapfile_no_mgrisae

while read proteome; do
    sbatch --job-name=$proteome.pannzer --mem=4000M -n 1 --export=ALL,proteome=$proteome \
        -p savio4_htc --qos=minium_htc4_normal --account=co_minium \
        /global/home/users/pierrj/git/slurm/pannzer_single_core.slurm
done < $MAPFILE


MAPFILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/proteomes_mapfile_no_mgrisae
while read proteome; do
    cp GO.${proteome}.out GO.${proteome}.out.copy
done < $MAPFILE

while read proteome; do
    wc -l GO.${proteome}.out
done < $MAPFILE

while read proteome; do
    protein_count=$(grep gene GO.${proteome}.out | awk '{print $1}' | sort | uniq | wc -l)
    b=7000
# grep gene GO.${proteome}.out | awk '{print $1}' | sort | uniq | wc -l
    if (( protein_count < b )); then
        echo $proteome
    fi
done < $MAPFILE

# srun -p savio4_htc -n 1 --qos=minium_htc4_normal --account=co_minium -t 60 --pty bash
# sbatch -p savio4_htc --qos=minium_htc4_normal --account=co_minium -n 56

# sbatch --mem=510000M -n 224 -p savio4_htc --qos=minium_htc4_normal --account=co_minium slurm/expression_guy11.slurm


## to make a single file for go terms in the random forest model

if [ -f "GO.all.out" ]; then
    rm "GO.all.out"
fi

while read proteome; do
    grep gene GO.${proteome}.out | awk '{print $1}' | sort | uniq >> GO.all.out
done < $MAPFILE






## AND FOR WHEAT BLAST TOOO

cd /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all/random_forest/go

PROTEOMES_PATH=/global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all/all_proteomes_corrected/
MAPFILE=/global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all/proteomes_mapfile_no_mgrisae

while read proteome; do
    sbatch --job-name=$proteome.pannzer --mem=4000M -n 1 --export=ALL,PROTEOMES_PATH=$PROTEOMES_PATH,proteome=$proteome \
        -p savio4_htc --qos=minium_htc4_normal --account=co_minium \
        /global/home/users/pierrj/git/slurm/pannzer_single_core.slurm
done < $MAPFILE