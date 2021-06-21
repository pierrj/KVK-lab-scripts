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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav_validation

GENE_GFF=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.gff3


python /global/home/users/pierrj/git/python/pav_validation_get_all_gene_deletions.py guy11_vs_others ${GENE_GFF} gffs output_sizes all_gene_deletions_output

source activate mummer
module load gnuplot
module load imagemagick

cd all_gene_deletions_output

REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/guy11_genome_baoetal2017.fasta
GENOME_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes
SUBDIR_OUTPUT=output_jpgs
PERCENT_ZEROES_FILTER=1
ALL_JPGS_DIR=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav_validation/jpgs_small_deletions
ALL_COORDS_DIR=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/pav_validation/coords_small_deletions
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
    /global/home/users/pierrj/git/bash/auto_mummer_plot.sh -r ${REFERENCE} -g ${GENOME_DIR} -s ${SUBDIR_OUTPUT} -p ${PERCENT_ZEROES_FILTER}
    cd ..
done

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


python /global/home/users/pierrj/git/python/pav_validation_write_confirmed_deletions.py ${ALL_COORDS_DIR} ${ALL_GENELOCS_DIR} all_small_deletions.bed

awk -v OFS='\t' '{print $1, $2, $3}' all_small_deletions.bed | sort -k1,1 -k2,2n | uniq > uniq.all_small_deletions.bed