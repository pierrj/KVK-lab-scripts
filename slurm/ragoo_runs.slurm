#!/bin/bash
#SBATCH --job-name=ragoo_runs
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
cd /global/scratch/users/pierrj/ragoo_runs
export PYTHONPATH=$PYTHONPATH:/global/scratch/users/pierrj/scripts/ragoo_python_install//lib/python3.7/site-packages/
source activate ragoo

while read genome
do
cd ${genome}
rm -r ragoo_output
cp /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/results/genes_on_eccs/G3_all.common.genes.gff3 .
python /global/scratch/users/pierrj/scripts/RaGOO/ragoo.py guy11_genome_baoetal2017.fasta moryzae_${genome}_ref.fasta -m /global/scratch/users/pierrj/software/minimap2/minimap2-2.17_x64-linux/minimap2 -b -gff G3_all.common.genes.gff3
ls ragoo_output/orderings/* > orderings.fofn
samtools faidx ragoo_output/chimera_break/guy11_genome_baoetal2017.intra.chimera.broken.fa
python /global/scratch/users/pierrj/scripts/RaGOO/lift_over.py ragoo_output/chimera_break/G3_all.common.genes.intra.chimera_broken.gff orderings.fofn ragoo_output/chimera_break/guy11_genome_baoetal2017.intra.chimera.broken.fa.fai > genes.ragoo.gff
cp genes.ragoo.gff ${genome}.liftover.gff
cp ragoo_output/ragoo.fasta ${genome}.ragoo.fasta
cd ..
done < mapfile_nobreak

while read genome
do
cd ${genome}_nobreak
rm -r ragoo_output
cp /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/results/genes_on_eccs/G3_all.common.genes.gff3 .
python /global/scratch/users/pierrj/scripts/RaGOO/ragoo.py guy11_genome_baoetal2017.fasta moryzae_${genome}_ref.fasta -m /global/scratch/users/pierrj/software/minimap2/minimap2-2.17_x64-linux/minimap2
ls ragoo_output/orderings/* > orderings.fofn
samtools faidx guy11_genome_baoetal2017.fasta
python /global/scratch/users/pierrj/scripts/RaGOO/lift_over.py G3_all.common.genes.gff3 orderings.fofn guy11_genome_baoetal2017.fasta.fai > genes.ragoo.gff
cp genes.ragoo.gff ${genome}_nobreak.liftover.gff
cp ragoo_output/ragoo.fasta ${genome}_nobreak.ragoo.fasta
cd ..
done < mapfile_nobreak