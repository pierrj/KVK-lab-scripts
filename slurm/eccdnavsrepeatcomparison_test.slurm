#!/bin/bash
#SBATCH --job-name=eccdnavsrepeatcomparison_test
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/repeat_comparison_wrapper/

GENOME_FASTA=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
SAMPLE=G3
THREADS=${SLURM_NTASKS}
REPEAT_FILE=/global/scratch/users/pierrj/references/te_annotations/moryzae/moryzae_non_ltr_repeats.bed
SAMPLE_MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/rnaseq_comparison_wrapper/anyoverlap_normalization/G3.sample_mapfile



ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/repeats_v_eccdnas/large_eccdnas/G3.large_eccdnas_mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do
cd ${sample}
    realpath large_eccdnas_splitreads >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/repeats_v_eccdnas/large_eccdnas/

/global/home/users/pierrj/git/bash/ecc_dna_vs_repetitive.sh -g ${GENOME_FASTA} \
    -s ${SAMPLE} -t ${THREADS} -f ${REPEAT_FILE} \
    -e ${ECCDNA_MAPFILE} -m ${SAMPLE_MAPFILE}



ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/repeats_v_eccdnas/micro_dnas/G3.micro_dnas_mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do
cd ${sample}
    realpath micro_dnas_splitreads >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/repeats_v_eccdnas/micro_dnas/

/global/home/users/pierrj/git/bash/ecc_dna_vs_repetitive.sh -g ${GENOME_FASTA} \
    -s ${SAMPLE} -t ${THREADS} -f ${REPEAT_FILE} \
    -e ${ECCDNA_MAPFILE} -m ${SAMPLE_MAPFILE}