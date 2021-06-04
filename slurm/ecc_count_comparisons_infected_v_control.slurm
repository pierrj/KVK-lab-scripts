#!/bin/bash
#SBATCH --job-name=ecc_count_comparisons_infected_v_control
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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/comparative/infected_v_control

if [ -f "ecc_count_comparison.txt" ]; then
    rm ecc_count_comparison.txt
fi

while read -r dir genome_file gff_file copia_file gypsy_file mapfile; do
    cd $dir
    while read sample;
    do
        cd ${sample}
        echo $sample 1>&2
        awk '{if ($3 == "gene") print $0}' ${gff_file} > tmp.genes.gff
        awk ' $5 != "lowq" ' ${sample}.ecc_caller_out.details.txt > ${sample}.ecc_caller_out.details.nolowq.txt
        bedtools intersect -u -f 0.5 -a ${sample}.ecc_caller_out.details.nolowq.txt -b tmp.genes.gff > tmp.genic_eccs
        bedtools intersect -v -f 0.5 -a ${sample}.ecc_caller_out.details.nolowq.txt -b tmp.genes.gff > tmp.noncoding_eccs
        bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b ${copia_file} | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) >= 0.9) {print $0}}' | sort -k1,1 -k2,2n > tmp.copia_eccs
        bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b ${gypsy_file} | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) >= 0.9) {print $0}}' | sort -k1,1 -k2,2n > tmp.gypsy_eccs
        bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b ${copia_file} | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $0}}' | \
            bedtools intersect -wao -a - -b ${gypsy_file} | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $0}}'| sort -k1,1 -k2,2n | awk '$3-$2 <= 400' > tmp.micro_dnas
        bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b ${copia_file} | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $0}}' | \
            bedtools intersect -wao -a - -b ${gypsy_file} | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $0}}'| sort -k1,1 -k2,2n | awk '$3-$2 > 400' > tmp.large_eccdnas
        ecc_count_genic=$(wc -l tmp.genic_eccs | awk '{print $1}')
        ecc_count_noncoding=$(wc -l tmp.noncoding_eccs | awk '{print $1}')
        ecc_count_copia=$(wc -l tmp.copia_eccs | awk '{print $1}')
        ecc_count_gypsy=$(wc -l tmp.gypsy_eccs | awk '{print $1}')
        ecc_count_micro=$(wc -l tmp.micro_dnas | awk '{print $1}')
        ecc_count_large=$(wc -l tmp.large_eccdnas | awk '{print $1}')
        ecc_count_all=$(wc -l ${sample}.ecc_caller_out.details.nolowq.txt | awk '{print $1}')
        read_count=$(samtools view -c no_secondary.filtered.sorted.${sample}.bam | awk '{print $1}') ## mapped reads
        genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum;}' ${genome_file} ) ## basepairs
        coding_bp=$(awk '{SUM+=$5-$4}END{print SUM}' tmp.genes.gff)
        copia_bp=$(awk '{SUM+=$3-$2}END{print SUM}' ${copia_file})
        gypsy_bp=$(awk '{SUM+=$3-$2}END{print SUM}' ${gypsy_file})
        echo -e $sample'\t'$ecc_count_genic'\t'$ecc_count_noncoding'\t'$ecc_count_copia\
            '\t'$ecc_count_gypsy'\t'$ecc_count_micro'\t'$ecc_count_large'\t'$ecc_count_all\
            '\t'$read_count'\t'$genome_size\
            '\t'$coding_bp'\t'$copia_bp'\t'$gypsy_bp >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/comparative/infected_v_control/ecc_count_comparison.txt
        rm tmp.genes.gff
        rm tmp.genic_eccs
        rm tmp.noncoding_eccs
        rm tmp.copia_eccs
        rm tmp.gypsy_eccs
        rm tmp.micro_dnas
        rm tmp.large_eccdnas
        cd ..
    done < $mapfile
done < mapfile_ecc_count