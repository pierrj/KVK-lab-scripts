#!/bin/bash
#SBATCH --job-name=acs_finding_final
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/acs

GENOME_FASTA=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
GENOME_CHROMSIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes
LTR_TE_LOCATIONS=/global/scratch/pierrj/references/te_annotations/moryzae/moryzae_copiaandgypsy_locs.bed

module load meme

if [ -f "large_eccdnas_all" ]; then
    rm large_eccdnas_all
fi

while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/large_eccdnas >> large_eccdnas_all; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
sort -k1,1 -k2,2n large_eccdnas_all | uniq > large_eccdnas_all_uniq

if [ -d fimo_out ]; then
    rm -r fimo_out
fi

bedtools getfasta -fi ${GENOME_FASTA} -bed large_eccdnas_all_uniq > large_eccdnas_all_uniq.fasta
fimo --verbosity 1 -o fimo_out acs_meme.txt large_eccdnas_all_uniq.fasta
## pick score cutoff
awk '$1!="#" && $6>=17' fimo_out/fimo.txt | awk -v OFS='\t' '{ split($2,a,":"); split(a[2],b,"-"); print a[1], b[1], b[2] }' | sort | uniq > large_eccdnas_w_acs
awk 'NR==FNR{a[$0]++;next}!a[$0]' large_eccdnas_w_acs large_eccdnas_all_uniq > large_eccdnas_no_acs

if [ -f "large_eccdnas_splitreads_all" ]; then
    rm large_eccdnas_splitreads_all
fi

while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/large_eccdnas_splitreads >> large_eccdnas_splitreads_all; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
sed 's/[[:space:]]*$//' large_eccdnas_splitreads_all | sort -k1,1 -k2,2n  | uniq -c | awk -v OFS='\t' '{print $2, $3, $4, $1}' > large_eccdnas_splitreads_all_counted

awk 'NR==FNR{a[$0]++;next}a[$1"\t"$2"\t"$3]' large_eccdnas_w_acs large_eccdnas_splitreads_all_counted > large_eccdnas_w_acs_counted
awk 'NR==FNR{a[$0]++;next}a[$1"\t"$2"\t"$3]' large_eccdnas_no_acs large_eccdnas_splitreads_all_counted > large_eccdnas_no_acs_counted

if [ -f permute_large_eccs.txt ]; then
    rm permute_large_eccs.txt
fi
for i in {1..100}; do 
    if [ -d shuffled_fimo ]; then
        rm -r shuffled_fimo
    fi
    bedtools shuffle -i large_eccdnas_all_uniq -g ${GENOME_CHROMSIZES} -excl ${LTR_TE_LOCATIONS} > shuffled
    bedtools getfasta -fi ${GENOME_FASTA} -bed shuffled > shuffled.fasta
    fimo --verbosity 1 -o shuffled_fimo acs_meme.txt shuffled.fasta
    ## pick score cutoff
    awk '$6>=17' shuffled_fimo/fimo.txt | awk '{print $2}' | sort | uniq | wc -l >> permute_large_eccs.txt
done







if [ -f "micro_dnas_all" ]; then
    rm micro_dnas_all
fi

while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/micro_dnas >> micro_dnas_all; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
sort -k1,1 -k2,2n micro_dnas_all | uniq > micro_dnas_all_uniq

if [ -d fimo_out ]; then
    rm -r fimo_out
fi

bedtools getfasta -fi ${GENOME_FASTA} -bed micro_dnas_all_uniq > micro_dnas_all_uniq.fasta
fimo --verbosity 1 -o fimo_out acs_meme.txt micro_dnas_all_uniq.fasta
## pick score cutoff
awk '$1!="#" && $6>=17' fimo_out/fimo.txt | awk -v OFS='\t' '{ split($2,a,":"); split(a[2],b,"-"); print a[1], b[1], b[2] }' | sort | uniq > micro_dnas_w_acs
awk 'NR==FNR{a[$0]++;next}!a[$0]' micro_dnas_w_acs micro_dnas_all_uniq > micro_dnas_no_acs

if [ -f "micro_dnas_splitreads_all" ]; then
    rm micro_dnas_splitreads_all
fi

while read sample; do cat /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/micro_dnas_splitreads >> micro_dnas_splitreads_all; done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
sed 's/[[:space:]]*$//' micro_dnas_splitreads_all | sort -k1,1 -k2,2n  | uniq -c | awk -v OFS='\t' '{print $2, $3, $4, $1}' > micro_dnas_splitreads_all_counted

awk 'NR==FNR{a[$0]++;next}a[$1"\t"$2"\t"$3]' micro_dnas_w_acs micro_dnas_splitreads_all_counted > micro_dnas_w_acs_counted
awk 'NR==FNR{a[$0]++;next}a[$1"\t"$2"\t"$3]' micro_dnas_no_acs micro_dnas_splitreads_all_counted > micro_dnas_no_acs_counted

if [ -f permute_micro_dnas.txt ]; then
    rm permute_micro_dnas.txt
fi
for i in {1..100}; do 
    if [ -d shuffled_fimo ]; then
        rm -r shuffled_fimo
    fi
    bedtools shuffle -i micro_dnas_all_uniq -g ${GENOME_CHROMSIZES} -excl ${LTR_TE_LOCATIONS} > shuffled
    bedtools getfasta -fi ${GENOME_FASTA} -bed shuffled > shuffled.fasta
    fimo --verbosity 1 -o shuffled_fimo acs_meme.txt shuffled.fasta
    ## pick score cutoff
    awk '$6>=17' shuffled_fimo/fimo.txt | awk '{print $2}' | sort | uniq | wc -l >> permute_micro_dnas.txt
done