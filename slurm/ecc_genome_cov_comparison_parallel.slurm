#!/bin/bash
#SBATCH --job-name=ecc_genome_cov_comparison_parallel
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

cd $ECC_DIR

cd ${sample}

GENOME_FILE=$(realpath ${GENOME_FILE}.filtered)
GENOME_CHROMSIZES=$(realpath ${GENOME_FILE}.chromsizes)

awk ' $5 != "lowq" ' ${sample}.ecc_caller_out.details.txt > ${sample}.ecc_caller_out.details.nolowq.txt
## base pairs
genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum;}' ${GENOME_FILE} )
## coverage base pairs
covered_bps=$(bedtools genomecov -d -i ${sample}.ecc_caller_out.details.nolowq.txt -g ${GENOME_FILE}.chromsizes | awk '$3 == 0' | wc -l)
echo -e $sample'\t'$genome_size'\t'$covered_bps'\t'$count >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/comparative/genome_coverage/ecc_coverage_comparison.txt