#!/bin/bash
#SBATCH --job-name=compare_q_option
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/compare_q_option/noqoption

FILTERED_BAMFILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/filtered.sorted.G3_1A.bam

/global/home/users/pierrj/git/bash/get_splitreads_for_LTR.sh -s G3_1A -b ${FILTERED_BAMFILE}

samtools view SRs.no_orientation.allmapq.G3_1A.bam | awk -v OFS='\t' '{print $1, $5}' | sort -k1,1 > SRs_and_qscore_sorted


awk -v OFS='\t' '{
    prev=$0; f1=$1; f2=$2
    getline 
    if ($2 != 0 && f2 != 0) {
        print "no_mapq0"
    }
    else if ($2 != 0 && f2 == 0) {
        print "one_mapq0"
    }
    else if ($2 == 0 && f2 != 0) {
        print "one_mapq0"
    }
    else if ($2 == 0 && f2 == 0) {
        print "both_mapq0"
    }
}'  SRs_and_qscore_sorted | sort | uniq -c > count


cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/compare_q_option/yesqoption

FILTERED_BAMFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_calling_with_mapq0/G3_1A_multimapping_qoption/filtered.sorted.G3_1A.bam

/global/home/users/pierrj/git/bash/get_splitreads_for_LTR.sh -s G3_1A -b ${FILTERED_BAMFILE}


samtools view SRs.no_orientation.allmapq.G3_1A.bam | awk -v OFS='\t' '{print $1, $5}' | sort -k1,1 > SRs_and_qscore_sorted

awk -v OFS='\t' '{
    prev=$0; f1=$1; f2=$2
    getline 
    if ($2 != 0 && f2 != 0) {
        print "no_mapq0"
    }
    else if ($2 != 0 && f2 == 0) {
        print "one_mapq0"
    }
    else if ($2 == 0 && f2 != 0) {
        print "one_mapq0"
    }
    else if ($2 == 0 && f2 == 0) {
        print "both_mapq0"
    }
}'  SRs_and_qscore_sorted | sort | uniq -c > count