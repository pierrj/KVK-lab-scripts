#!/bin/bash
#SBATCH --job-name=ecc_sv_finding_across_family
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_sv_finding_in_sordariomycetidae/


GENOMES=/global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_sv_finding_in_sordariomycetidae/
MORYZAE_GENOME=/global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_sv_finding_in_sordariomycetidae/Magor1_AssemblyScaffolds_Repeatmasked.fasta
PERCENT_ZEROES_FILTER=1.0

ls -1 | grep -v mapfile_taxa | awk '{print substr($1,0,length($1)-37)}' > mapfile_taxa

if [ -f "jobqueue" ]; then
    rm "jobqueue"
fi

while read taxon;
do
    TAXON_GENOME=${GENOMES}/${taxon}_AssemblyScaffolds_Repeatmasked.fasta
    echo "/global/home/users/pierrj/git/bash/ecc_sv_finder.sh -r ${MORYZAE_GENOME} -q ${TAXON_GENOME} -t 1 -n moryzae_v_${taxon} -p ${PERCENT_ZEROES_FILTER} -o ${GENOMES}/moryzae_v_${taxon}_jpg -m ${GENOMES}/moryzae_v_${taxon}_tmp" >> jobqueue
    echo "/global/home/users/pierrj/git/bash/ecc_sv_finder.sh -r ${TAXON_GENOME} -q ${MORYZAE_GENOME} -t 1 -n ${taxon}_v_moryzae -p ${PERCENT_ZEROES_FILTER} -o ${GENOMES}/${taxon}_v_moryzae_jpg -m ${GENOMES}/${taxon}_v_moryzae_tmp" >> jobqueue
done < mapfile_taxa

parallel -j ${SLURM_NTASKS} < jobqueue