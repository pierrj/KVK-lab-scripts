#!/bin/bash
#SBATCH --job-name=ecc_count_comparisons_codingnoncoding
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
cd /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/

while read -r dir ecc_count_output genome_file gff_file mapfile; do
    cd $dir
    echo $ecc_count_output 1>&2
    ecc_count_output_genic=ecc_count_${ecc_count_output}_genic.txt
    ecc_count_output_noncoding=ecc_count_${ecc_count_output}_noncoding.txt
    ecc_count_output_all=ecc_count_${ecc_count_output}_all.txt
    if [ -f "ecc_count_${ecc_count_output}_genic.txt" ]; then
        rm ecc_count_${ecc_count_output}_genic.txt
    fi
    if [ -f "ecc_count_${ecc_count_output}_noncoding.txt" ]; then
        rm ecc_count_${ecc_count_output}_noncoding.txt
    fi
    if [ -f "ecc_count_${ecc_count_output}_all.txt" ]; then
        rm ecc_count_${ecc_count_output}_all.txt
    fi   
    while read sample;
    do
        cd ${sample}
        echo $sample 1>&2
        awk '{if ($3 == "gene") print $0}' ${gff_file} > tmp.genes.gff
        bedtools intersect -u -f 0.5 -a ecccaller_output.${sample}.renamed.details.tsv -b tmp.genes.gff > tmp.genic.ecccaller_output.${sample}.renamed.details.tsv
        bedtools intersect -v -f 0.5 -a ecccaller_output.${sample}.renamed.details.tsv -b tmp.genes.gff > tmp.noncoding.ecccaller_output.${sample}.renamed.details.tsv
        ecc_count_genic=$(awk ' $6 != "lowq" ' tmp.genic.ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
        ecc_count_noncoding=$(awk ' $6 != "lowq" ' tmp.noncoding.ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
        ecc_count_all=$(awk ' $6 != "lowq" ' ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
        read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
        genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${genome_file} ) ## megabasepairs
        total_genic=$(awk -v ECC="$ecc_count_genic" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
        total_noncoding=$(awk -v ECC="$ecc_count_noncoding" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
        total_all=$(awk -v ECC="$ecc_count_all" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
        echo -e $sample'\t'$total_genic >> ../ecc_count_output_genic
        echo -e $sample'\t'$total_noncoding >> ../ecc_count_output_noncoding
        echo -e $sample'\t'$total_all >> ../ecc_count_output_all
        rm tmp.genes.gff
        rm tmp.genic.ecccaller_output.${sample}.renamed.details.tsv
        rm tmp.noncoding.ecccaller_output.${sample}.renamed.details.tsv
        cd ..
    done < $mapfile
    cp ${ecc_count_output} /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/
done < mapfile_ecc_count