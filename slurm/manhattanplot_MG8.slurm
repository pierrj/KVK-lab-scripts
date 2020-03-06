#!/bin/bash
#SBATCH --job-name=manhattanplot_G3_1A_MG8_color
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/SMRTcell_all/G3_1A
R CMD BATCH manhattanplot_MG8.R
