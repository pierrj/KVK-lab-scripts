#!/bin/bash
#SBATCH --job-name=deletion_mummer_validation
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

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/deletion_mummer_plots/guy11_v_all_redo_with_qc/gene_deletion_output

source activate mummer
module load gnuplot
module load imagemagick

REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/guy11_genome_baoetal2017.fasta
GENOME_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes
SUBDIR_OUTPUT=output_jpgs
PERCENT_ZEROES_FILTER=1

for gene_dir in */ ; do
    gene=$(basename ${gene_dir})
    cd ${gene}
    /global/home/users/pierrj/git/bash/auto_mummer_plot.sh -r ${REFERENCE} -g ${GENOME_DIR} -s ${SUBDIR_OUTPUT} -p ${PERCENT_ZEROES_FILTER}
    cd ..
done

# for d in */ ; do
#     cd $d
#         d_basename=$(basename $d)
#             cd output_jpgs
#             for f in *; do
#                 cp ${f} /global/scratch/users/pierrj/eccDNA/pipeline_tests/deletion_mummer_plots/guy11_v_all/all_jpgs/${d_basename}_${f}
#             done
#         cd ..
#     cd ..
# done

for gene_dir in */ ; do
    gene=$(basename ${gene_dir})
    cd ${gene}
    while read genome; do
        cd guy11_v_${genome}
            cd tmp
                show-coords 0.delta > 0.coords
                tail -n+6 0.coords | awk -v OFS='\t' '{print $1, $2, $4, $5, $12, $13}' > ${gene}_${genome}.processed_output.coords
                python /global/home/users/pierrj/git/python/process_candidate_deletions.py ${gene}_${genome}.processed_output.coords ../guy11_geneloc.bed ../${gene}_${genome}.deletion.bed 40 ${gene} ${genome}
                if [ -f "../${gene}_${genome}.deletion.bed" ]; then
                    cp ../${gene}_${genome}.deletion.bed /global/scratch/users/pierrj/eccDNA/pipeline_tests/deletion_mummer_plots/guy11_v_all/all_deletions/${gene}_${genome}.deletion.bed
                fi
            cd ..
        cd ..
    done < genome_mapfile
    cd ..
done