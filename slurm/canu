#!/bin/bash
#SBATCH --job-name=canu
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/SMRTcell_all/G3_1A/canu/meta
canu -assemble -p metagenome_test -pacbio-corrected ../demux.fastq genomeSize=43m useGrid=false corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=200
