#!/bin/bash
#SBATCH --job-name=minichromosome_exonerate
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

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/minichromosomes

module load exonerate

isolate=FR13

exonerate --model protein2genome guy11_fungap_out_prot.faa ${isolate}_no_minichrom.fasta --showtargetgff yes --bestn 1 --showalignment no --verbose 0 --showvulgar no > ${isolate}_genes_core.gff

awk '$3=="gene"' ${isolate}_genes_core.gff | awk '$19>70' | awk '{print $13}' > ${isolate}_genes_on_core

isolate=CD156

exonerate --model protein2genome guy11_fungap_out_prot.faa ${isolate}_no_minichrom.fasta --showtargetgff yes --bestn 1 --showalignment no --verbose 0 --showvulgar no > ${isolate}_genes_core.gff

awk '$3=="gene"' ${isolate}_genes_core.gff | awk '$19>70' | awk '{print $13}' > ${isolate}_genes_on_core

isolate=US71

exonerate --model protein2genome guy11_fungap_out_prot.faa ${isolate}_no_minichrom.fasta --showtargetgff yes --bestn 1 --showalignment no --verbose 0 --showvulgar no > ${isolate}_genes_core.gff

awk '$3=="gene"' ${isolate}_genes_core.gff | awk '$19>70' | awk '{print $13}' > ${isolate}_genes_on_core