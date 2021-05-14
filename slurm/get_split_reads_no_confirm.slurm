#!/bin/bash
#SBATCH --job-name=get_split_reads_no_confirm
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

MAPFILE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames
SAMPLE=G3_1A
THREADS=20
FILTERED_BAMFILE=uniq.filtered.sorted.G3_1A.bam


cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A

samtools view -f 81 -F 4 ${FILTERED_BAMFILE} > tmp.reverseread1.${SAMPLE}.sam
splitread_file="reverseread1.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -f 145 -F 4 ${FILTERED_BAMFILE} > tmp.reverseread2.${SAMPLE}.sam
splitread_file="reverseread2.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -f 65 -F 20 ${FILTERED_BAMFILE} > tmp.forwardread1.${SAMPLE}.sam
splitread_file="forwardread1.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -f 129 -F 20 ${FILTERED_BAMFILE} > tmp.forwardread2.${SAMPLE}.sam
splitread_file="forwardread2.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

# I did some merging of reads in the past but this ended up being detrimental
# these next two chunks should be removed eventually
# currently both of these should contain 0 reads
samtools view -f 16 -F 5 ${FILTERED_BAMFILE} > tmp.reversemerged.${SAMPLE}.sam
splitread_file="reversemerged.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -F 21 ${FILTERED_BAMFILE} > tmp.forwardmerged.${SAMPLE}.sam
splitread_file="forwardmerged.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

# putting them all back together
cat tmp.oriented.samechromosome.exactlytwice.qualityfiltered.reverseread1.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.reverseread2.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.forwardread1.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.forwardread2.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.reversemerged.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.forwardmerged.${SAMPLE}.sam > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam

# verify that the match lengths on either side match up so that only proper split reads are looked on
python /global/home/users/pierrj/git/python/filter_for_match_lengths.py tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam tmp.match_length_filtered.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam

# converting to bed file
samtools view -b -h <(cat <(samtools view -H ${FILTERED_BAMFILE}) tmp.match_length_filtered.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam) > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.bam
bedtools bamtobed -i tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.bam | sort -k4,4 -k2,2n > splitreads.${SAMPLE}.bed

# merging split read halves into single, putative eccDNA forming regions to be confirmed or rejected
awk -v OFS='\t' '{
    prev=$0; f2=$2; f4=$4
    getline 
    if ($4 == f4 && f2 < $2) {
        print $1, f2, $3, $4
    }
}' splitreads.${SAMPLE}.bed > merged.splitreads.${SAMPLE}.bed

# length filter because we don't expect eccDNAs to be that big
# could be tweaked potentially but this gets rid of very few split reads
# 50k is based off the approximate size of eccDNAs that should be coming out of the column NEEDS TO BE VERIFIED
awk -v OFS='\t' '$3-$2<50000' merged.splitreads.${SAMPLE}.bed > lengthfiltered.merged.splitreads.${SAMPLE}.bed