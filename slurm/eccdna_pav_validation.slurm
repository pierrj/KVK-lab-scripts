#!/bin/bash
#SBATCH --job-name=eccdna_pav_validation
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav_validation

# gene annotation file
GENE_GFF=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.gff3

# this script looks for genes that are found in guy11 genomes that are deleted in other moryzae genomes
# specifically, these must be single gene deletions, where the other moryzae genomes have kept the neighboring genes
# it outputs the location of putative deletions in other genomes
python /global/home/users/pierrj/git/python/pav_validation_get_all_gene_deletions.py guy11_vs_others ${GENE_GFF} gffs output_sizes all_gene_deletions_output

source activate mummer
module load gnuplot
module load imagemagick

cd all_gene_deletions_output

REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/guy11_genome_baoetal2017.fasta
# directory containing fasta files for all other moryzae genomes
GENOME_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes
SUBDIR_OUTPUT=output_jpgs
# for auto_mummer_plot.sh, this means no alignments are filtered out bc they are incomplete
PERCENT_ZEROES_FILTER=1
# directory where jpegs of mummer alignments are stored
ALL_JPGS_DIR=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav_validation/jpgs_small_deletions
# directory where nucmer alignment files are stored
ALL_COORDS_DIR=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav_validation/coords_small_deletions
# directory where the location of deletions in genomes are
ALL_GENELOCS_DIR=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav_validation/genelocs_small_deletions

if [ -d "${ALL_JPGS_DIR}" ]; then
    rm -r ${ALL_JPGS_DIR}
fi

if [ -d "${ALL_COORDS_DIR}" ]; then
    rm -r ${ALL_COORDS_DIR}
fi

if [ -d "${ALL_GENELOCS_DIR}" ]; then
    rm -r ${ALL_GENELOCS_DIR}
fi

mkdir ${ALL_JPGS_DIR}
mkdir ${ALL_COORDS_DIR}
mkdir ${ALL_GENELOCS_DIR}

for gene_dir in */ ; do
    gene=$(basename ${gene_dir})
    cd ${gene}
    # nucmer align locations of putative deletions to validate them
    /global/home/users/pierrj/git/bash/auto_mummer_plot.sh -r ${REFERENCE} -g ${GENOME_DIR} -s ${SUBDIR_OUTPUT} -p ${PERCENT_ZEROES_FILTER}
    cd ..
done

# move output jpgs from auto_mummer_plot.sh
for d in */ ; do
    cd $d
        d_basename=$(basename $d)
            cd output_jpgs
            for f in *; do
                cp ${f} ${ALL_JPGS_DIR}/${d_basename}_${f}
            done
        cd ..
    cd ..
done

# process outputs from auto_mummer_plot.sh for final validation
for gene_dir in */ ; do
    gene=$(basename ${gene_dir})
    cd ${gene}
    while read genome; do
        cd guy11_v_${genome}
            cd tmp
                show-coords 0.delta > 0.coords
                tail -n+6 0.coords | awk -v OFS='\t' '{print $1, $2, $4, $5, $12, $13}' > ${gene}_${genome}.processed_output.coords
                cp ${gene}_${genome}.processed_output.coords ${ALL_COORDS_DIR}/${gene}_${genome}.processed_output.coords
                cp ../guy11_geneloc.bed ${ALL_GENELOCS_DIR}/${gene}_${genome}.gene_loc
            cd ..
        cd ..
    done < genome_mapfile
    cd ..
done

# finally, validate locations of deletions by parsing nucmer outputs and making sure that a deletion is actually there in the genome and not some sort of misannotation
python /global/home/users/pierrj/git/python/pav_validation_write_confirmed_deletions.py ${ALL_COORDS_DIR} ${ALL_GENELOCS_DIR} all_small_deletions.bed

awk -v OFS='\t' '{print $1, $2, $3}' all_small_deletions.bed | sort -k1,1 -k2,2n | uniq > uniq.all_small_deletions.bed