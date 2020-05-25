#!/bin/bash
#SBATCH --job-name=sra_fastqdump
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
outdir=rnaseq/moryzae/2012_soanes
accession=SRR298847
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${accession} -O /global/scratch/users/pierrj/${outdir}
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e 24 -O /global/scratch/users/pierrj/${outdir} -t /global/scratch/users/pierrj/${outdir}/tmp /global/scratch/users/pierrj/${outdir}/${accession}.sra
outdir=rnaseq/moryzae/2012_soanes
accession=SRR298848
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${accession} -O /global/scratch/users/pierrj/${outdir}
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e 24 -O /global/scratch/users/pierrj/${outdir} -t /global/scratch/users/pierrj/${outdir}/tmp /global/scratch/users/pierrj/${outdir}/${accession}.sra
outdir=rnaseq/moryzae/2012_soanes
accession=SRR298849
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${accession} -O /global/scratch/users/pierrj/${outdir}
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e 24 -O /global/scratch/users/pierrj/${outdir} -t /global/scratch/users/pierrj/${outdir}/tmp /global/scratch/users/pierrj/${outdir}/${accession}.sra