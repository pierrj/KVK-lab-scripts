#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/blast
if [ $type == "n" ];
    then blastn -query /global/scratch/users/pierrj/references/${query} -db /global/scratch/users/pierrj/references/${db} -out ${out}.blastn
elif [ $type == "p" ];
    then blastp -query /global/scratch/users/pierrj/references/${query} -db /global/scratch/users/pierrj/references/${db} -out ${out}.blastp
elif [ $type == "x" ];
    then blastx -query /global/scratch/users/pierrj/references/${query} -db /global/scratch/users/pierrj/references/${db} -out ${out}.blastx
elif [ $type == "tn" ];
    then tblastn -query /global/scratch/users/pierrj/references/${query} -db /global/scratch/users/pierrj/references/${db} -out ${out}.tblastn
elif [ $type == "tx" ];
    then tblastx -query /global/scratch/users/pierrj/references/${query} -db /global/scratch/users/pierrj/references/${db} -out ${out}.tblastx
fi

### USAGE
### set your variables
### then use this command
### sbatch --export=type=$type,query=$query,db=$db,out=$out /global/home/users/pierrj/git/slurm/blast.slurm
### USAGE