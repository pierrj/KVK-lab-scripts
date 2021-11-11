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
#!/bin/bash
#SBATCH --job-name=ltr_splitreads_per_element_plot
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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/ltr_splitreads/per_element_plot

# bed file containing locations of te consensus sequences generated from wickersoft
TE_BEDFILE=/global/scratch/pierrj/references/te_annotations/moryzae/moryzae_high_conf_ltr_tes_nolengthfilter_individual_elements.bed
# mapfile of element names
ELEMENT_MAPFILE=/global/scratch/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/coverage_plots/transposon_mapfile
# mapfile of samples
MAPFILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
genome_file=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta

genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum;}' ${genome_file} )

# separate bed files by te
while read ELEMENT; do
    grep ${ELEMENT} ${TE_BEDFILE} > ${ELEMENT}.bed
done < ${ELEMENT_MAPFILE}

if [ -f "te_eccs_per_sample.txt" ]; then
    rm te_eccs_per_sample.txt
fi

while read sample; do
    if [ -f "tmp_line" ]; then
        rm tmp_line
    fi
    echo ${sample} >> tmp_line
    # count non-lowq eccs and split reads
    wc -l /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.ecc_caller_out.details.nolowq.txt >> tmp_line
    wc -l /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.ecc_caller_out.splitreads.bed >> tmp_line
    while read ELEMENT; do
        # intersect and count number of eccdns and split reads that are made up of 90% each element
        bedtools intersect -wao -a /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.ecc_caller_out.details.nolowq.txt -b ${ELEMENT}.bed | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) >= 0.9) {print $1, $2, $3}}' | sort -k1,1 -k2,2n > ${sample}.${ELEMENT}.eccs.bed
        sed 's/[[:space:]]*$//' /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.ecc_caller_out.splitreads.bed | awk 'NR==FNR{a[$0]++;next}a[$0]' ${sample}.${ELEMENT}.eccs.bed - > ${sample}.${ELEMENT}.splitreads.bed
        # count
        wc -l ${sample}.${ELEMENT}.eccs.bed | awk '{print $1}' >> tmp_line
        wc -l ${sample}.${ELEMENT}.splitreads.bed | awk '{print $1}' >> tmp_line
        # count bp of genome made up of ecc
        awk '{SUM+=$3-$2}END{print SUM}' ${ELEMENT}.bed >> tmp_line
    done < ${ELEMENT_MAPFILE}
    cat tmp_line | tr "\n" " " >> te_eccs_per_sample.txt
    echo >> te_eccs_per_sample.txt
done < ${MAPFILE}