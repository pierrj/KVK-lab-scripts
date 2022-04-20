
cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/comparative/genome_coverage

if [ -f "ecc_coverage_comparison.txt" ]; then
    rm ecc_coverage_comparison.txt
fi

let count=0

while read -r ORGANISM_NAME GENE_GFF SAMPLE_MAPFILE ECC_DIR GENOME_FILE COPIA_FILE GYPSY_FILE CONTIGNAMES_FILE; do
    cd $ECC_DIR
    while read sample; do
    let count=$count+1
    sbatch --job-name=$sample.genomecov --export=ECC_DIR=$ECC_DIR,sample=$sample,GENOME_FILE=$GENOME_FILE,count=$count\
         /global/home/users/pierrj/git/slurm/ecc_genome_cov_comparison_parallel.slurm
    done < $SAMPLE_MAPFILE
done < mapfile_ecc_count

cp ecc_coverage_comparison.txt ecc_coverage_comparison.txt.done

sort -k4,4n ecc_coverage_comparison.txt.done > ecc_coverage_comparison.txt.sorted

cut -f 1,2,3 ecc_coverage_comparison.txt.sorted > ecc_coverage_comparison.txt