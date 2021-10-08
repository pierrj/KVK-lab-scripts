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
THREADS=20


cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pacbio/all_false_positives


while read sample; do 
    cd ${sample}
        SAMPLE=$sample
        FILTERED_BAMFILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/uniq.filtered.sorted.${sample}.bam
        FILTERED_BAMFILE_QSORTED=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/multimapped.filtered.name_sorted.${sample}.bam

        samtools view -f 81 -F 4 ${FILTERED_BAMFILE} > tmp.reverseread1.${SAMPLE}.sam
        splitread_file="reverseread1.${SAMPLE}.sam"
        awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 9 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 9) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
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
        awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 9 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 9) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
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
        awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 9 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 9) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
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
        awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 9 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 9) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
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
        awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 9 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 9) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
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
        awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 9 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 9) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
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
        # python /global/home/users/pierrj/git/python/filter_for_match_lengths.py tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam tmp.match_length_filtered.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam

        cp tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam tmp.match_length_filtered.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam

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

        # get outward facing read pairs using sam flags
        # convert to bed file
        # fix names for filtering
        # filter to appearing only exactly twice, meaning that only complete read pairs are present
        samtools view ${FILTERED_BAMFILE} | awk '{ if (($2 == 81 || $2 == 83 || $2 == 145 || $2 == 147 ) && $9 > 0) print $0 ; else if (($2 == 97 || $2 == 99 || $2 == 161 || $2 == 163) && $9 <0) print $0}' | cat <(samtools view -H ${FILTERED_BAMFILE}) - | samtools view -b -h - > tmp.outwardfacing.${SAMPLE}.bam
        bedtools bamtobed -i tmp.outwardfacing.${SAMPLE}.bam | sort -k 4,4 > tmp.outwardfacing.${SAMPLE}.bed
        mv tmp.outwardfacing.${SAMPLE}.bed tmp.outwardfacing.${SAMPLE}.bed.old
        awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,substr($4, 1, length($4)-2),$5,$6}' tmp.outwardfacing.${SAMPLE}.bed.old > tmp.outwardfacing.${SAMPLE}.bed.old.trimmed
        awk 'NR==FNR{a[$4]++; next} a[$4]==2' tmp.outwardfacing.${SAMPLE}.bed.old.trimmed tmp.outwardfacing.${SAMPLE}.bed.old.trimmed > outwardfacing.${SAMPLE}.bed

        # change names of scaffolds using mapfiles for compatability with any genome
        chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
        for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > tmp.chrom_count
        paste tmp.chrom_count ${MAPFILE} > tmp.chrom_count_and_names
        awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names outwardfacing.${SAMPLE}.bed > outwardfacing.${SAMPLE}.renamed.bed
        awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names lengthfiltered.merged.splitreads.${SAMPLE}.bed > lengthfiltered.merged.splitreads.${SAMPLE}.renamed.bed

        # merge outward facing read pairs into single lines for confirming using python script
        sort -k4,4 -k2,2n outwardfacing.${SAMPLE}.renamed.bed > sorted.outwardfacing.${SAMPLE}.renamed.bed
        awk -v OFS='\t' '{
            prev=$0; f2=$2; f3=$3; f4=$4
            getline 
            if ($4 == f4 && f2 < $2 && f3 <$3) {
                print $1, f2, $3, f3, $2
            }
        }' sorted.outwardfacing.${SAMPLE}.renamed.bed > sorted.grouped.outwardfacing.${SAMPLE}.renamed.bed

        # use GNU parallel to speed things up
        # split into chunks first then each thread works on a chunk
        # parallel.confirmed are split reads confirmed by opposite facing read pairs
        split --number=l/${THREADS} --numeric-suffixes=1 --additional-suffix=.bed lengthfiltered.merged.splitreads.${SAMPLE}.renamed.bed lengthfiltered.merged.splitreads.${SAMPLE}.renamed.
        parallel -j ${THREADS} --link python /global/home/users/pierrj/git/python/ecc_caller_anygenome_confirmsrs_numpy_gnuparallel.py lengthfiltered.merged.splitreads.${SAMPLE}.renamed.{}.bed sorted.grouped.outwardfacing.${SAMPLE}.renamed.bed ${SAMPLE} ${chrom_count} {} ::: $(seq -w 1 ${THREADS})
        cat $(find . -maxdepth 1 -name "parallel.confirmed*" | xargs -r ls -1 | tr "\n" " ") > parallel.confirmed

        mv parallel.confirmed unique_parallel.confirmed

        rm parallel.confirmed*

        ## get length distribution of eccDNAs from uniquely mapped reads

        awk '{print $3-$2}' unique_parallel.confirmed > dsn.unique_parallel.confirmed

        ## primary only bc of the exactly twice filter
        # for multimapped read we want to broaden our scope of split reads
        # here we are including any split read that matches the required quality filter regardless of how they are oriented to each other
        samtools view -b -F 256 ${FILTERED_BAMFILE_QSORTED} > primary_only.${SAMPLE}.sorted.mergedandpe.bwamem.bam

        samtools view -f 65 -F 4 primary_only.${SAMPLE}.sorted.mergedandpe.bwamem.bam > tmp.primary_only.filtered.sorted.allmapq.mapped.1.${SAMPLE}.sam
        file='1'
        awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 9 ) print $0; else if (b !~ /[DMIHS]/ && int(b) > 9) print $0}' tmp.primary_only.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam
        awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam

        samtools view -f 129 -F 4 primary_only.${SAMPLE}.sorted.mergedandpe.bwamem.bam > tmp.primary_only.filtered.sorted.allmapq.mapped.2.${SAMPLE}.sam
        file='2'
        awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 9 ) print $0; else if (b !~ /[DMIHS]/ && int(b) > 9) print $0}' tmp.primary_only.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam
        awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam tmp.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam > tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.${file}.${SAMPLE}.sam

        # sort split reads by how many ends are multi-mapped
        awk -v OFS='\t' '{
            prev=$0; f1=$1 ; f5=$5
            getline 
            if ($1 == f1 && $5 != 0 && f5 != 0) {
                print f1 > "tmp.doubleunique.readnames.1"
            }
            else if ($1 == f1 && $5 != 0 && f5 == 0) {
                print f1 > "tmp.singleunique.readnames.1"
            }
            else if ($1 == f1 && $5 == 0 && f5 != 0) {
                print f1 > "tmp.singleunique.readnames.1"
            }
            else if ($1 == f1 && $5 == 0 && f5 == 0) {
                print f1 > "tmp.doublemapq0.readnames.1"
            }
        }'  tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.1.${SAMPLE}.sam

        awk -v OFS='\t' '{
            prev=$0; f1=$1 ; f5=$5
            getline 
            if ($1 == f1 && $5 != 0 && f5 != 0) {
                print f1 > "tmp.doubleunique.readnames.2"
            }
            else if ($1 == f1 && $5 != 0 && f5 == 0) {
                print f1 > "tmp.singleunique.readnames.2"
            }
            else if ($1 == f1 && $5 == 0 && f5 != 0) {
                print f1 > "tmp.singleunique.readnames.2"
            }
            else if ($1 == f1 && $5 == 0 && f5 == 0) {
                print f1 > "tmp.doublemapq0.readnames.2"
            }
        }'  tmp.exactlytwice.qualityfiltered.filtered.sorted.allmapq.mapped.2.${SAMPLE}.sam


        # split bamfiles with all reads and all alignments into read one and read 2
        samtools view -b -f 65 -F 4 ${FILTERED_BAMFILE_QSORTED} > tmp.filtered.sorted.allmapq.mapped.1.${SAMPLE}.bam
        samtools view -b -f 129 -F 4 ${FILTERED_BAMFILE_QSORTED} > tmp.filtered.sorted.allmapq.mapped.2.${SAMPLE}.bam

        # use picard filter sam reads to extract split reads identified above and all their possible secondary alignnents
        java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar FilterSamReads \
                INPUT=tmp.filtered.sorted.allmapq.mapped.1.${SAMPLE}.bam \
                OUTPUT=${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bam \
                READ_LIST_FILE=tmp.singleunique.readnames.1 \
                FILTER=includeReadList \
                SORT_ORDER=unsorted

        java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar FilterSamReads \
                INPUT=tmp.filtered.sorted.allmapq.mapped.2.${SAMPLE}.bam \
                OUTPUT=${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.2.bam \
                READ_LIST_FILE=tmp.singleunique.readnames.2 \
                FILTER=includeReadList \
                SORT_ORDER=unsorted

        java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar FilterSamReads \
                INPUT=tmp.filtered.sorted.allmapq.mapped.1.${SAMPLE}.bam \
                OUTPUT=${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bam \
                READ_LIST_FILE=tmp.doublemapq0.readnames.1 \
                FILTER=includeReadList \
                SORT_ORDER=unsorted

        java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/picard.jar FilterSamReads \
                INPUT=tmp.filtered.sorted.allmapq.mapped.2.${SAMPLE}.bam \
                OUTPUT=${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.2.bam \
                READ_LIST_FILE=tmp.doublemapq0.readnames.2 \
                FILTER=includeReadList \
                SORT_ORDER=unsorted

        # convert to bed files and merge
        bedtools bamtobed -cigar -i ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bed
        bedtools bamtobed -cigar -i ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.2.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.2.bed

        bedtools bamtobed -cigar -i ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bed
        bedtools bamtobed -cigar -i ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.2.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.2.bed

        cat ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.1.bed ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.2.bed > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.bed
        cat ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.1.bed ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.2.bed > ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.bed


        # split files into chunks for parallelization and fix chunks so that a read name isnt present in two chunks
        split --number=l/${THREADS} --numeric-suffixes=1 --additional-suffix=.bed ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.bed ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.chunk.

        for i in $(seq -w 1 1 $THREADS); do
            python /global/home/users/pierrj/git/python/split_chunk_fixer.py ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.chunk.${i}.bed ${i}
        done

        seq -w 1 1 $((THREADS-1)) > tmp.seq
        seq -w 2 1 ${THREADS} > tmp.seq_plusone
        paste tmp.seq tmp.seq_plusone > tmp.seqs

        while IFS=$'\t' read -r i next; do
            cat ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.chunk.${i}.bed split_line_fix.${next} > multimapped_splitreads.${i}.bed
        done < tmp.seqs

        cp ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.doublemapq0.chunk.${THREADS}.bed multimapped_splitreads.${THREADS}.bed

        # run through python script which calls putative eccdna forming regions from multimapping reads based off the length distribution from uniquely mapped eccDNAs
        parallel -j ${THREADS} --link python /global/home/users/pierrj/git/python/ecc_calling_mapq0.py dsn.unique_parallel.confirmed  multimapped_splitreads.{}.bed 50000 {} ::: $(seq -w 1 ${THREADS})

        python /global/home/users/pierrj/git/python/ecc_calling_mapq0_singleunique.py dsn.unique_parallel.confirmed ${SAMPLE}.sorted.mergedandpe.bwamem.multimapped_splitreads.singleunique.bed 50000

        # merge all chunks together then confirm as before
        cat mapq0_choices.* singleunique_choices > mapq0_single_unique_choices.bed


        # reget opposite facing reads as before but with mapq0s as well this time
        samtools view primary_only.${SAMPLE}.sorted.mergedandpe.bwamem.bam | awk '{ if (($2 == 81 || $2 == 83 || $2 == 145 || $2 == 147 ) && $9 > 0) print $0 ; else if (($2 == 97 || $2 == 99 || $2 == 161 || $2 == 163) && $9 <0) print $0}' | cat <(samtools view -H ${FILTERED_BAMFILE}) - | samtools view -b -h - > tmp.outwardfacing.${SAMPLE}.bam
        bedtools bamtobed -i tmp.outwardfacing.${SAMPLE}.bam | sort -k 4,4 > tmp.outwardfacing.${SAMPLE}.bed
        mv tmp.outwardfacing.${SAMPLE}.bed tmp.outwardfacing.${SAMPLE}.bed.old
        awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,substr($4, 1, length($4)-2),$5,$6}' tmp.outwardfacing.${SAMPLE}.bed.old > tmp.outwardfacing.${SAMPLE}.bed.old.trimmed
        awk 'NR==FNR{a[$4]++; next} a[$4]==2' tmp.outwardfacing.${SAMPLE}.bed.old.trimmed tmp.outwardfacing.${SAMPLE}.bed.old.trimmed > multi_mapping.outwardfacing.${SAMPLE}.bed

        awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names multi_mapping.outwardfacing.${SAMPLE}.bed > multi_mapping.outwardfacing.${SAMPLE}.renamed.bed
        awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names mapq0_single_unique_choices.bed > mapq0_single_unique_choices.renamed.bed


        # merge outward facing read pairs into single lines for confirming using python script
        sort -k4,4 -k2,2n multi_mapping.outwardfacing.${SAMPLE}.renamed.bed > sorted.multi_mapping.outwardfacing.${SAMPLE}.renamed.bed
        awk -v OFS='\t' '{
            prev=$0; f2=$2; f3=$3; f4=$4
            getline 
            if ($4 == f4 && f2 < $2 && f3 <$3) {
                print $1, f2, $3, f3, $2
            }
        }' sorted.multi_mapping.outwardfacing.${SAMPLE}.renamed.bed > sorted.grouped.multi_mapping.outwardfacing.${SAMPLE}.renamed.bed

        split --number=l/${THREADS} --numeric-suffixes=1 --additional-suffix=.bed mapq0_single_unique_choices.renamed.bed mapq0_single_unique_choices.renamed.
        parallel -j ${THREADS} --link python /global/home/users/pierrj/git/python/ecc_caller_anygenome_confirmsrs_numpy_gnuparallel.py mapq0_single_unique_choices.renamed.{}.bed sorted.grouped.multi_mapping.outwardfacing.${SAMPLE}.renamed.bed ${SAMPLE} ${chrom_count} {} ::: $(seq -w 1 ${THREADS})
        cat $(find . -maxdepth 1 -name "parallel.confirmed*" | xargs -r ls -1 | tr "\n" " ") > parallel.confirmed

        mv parallel.confirmed mapq0_parallel.confirmed

        cat unique_parallel.confirmed mapq0_parallel.confirmed > parallel.confirmed

        # convert scaffolds to 1 index from 0 index
        # rename scaffolds in parallel.confirmed
        paste ${MAPFILE} tmp.chrom_count > tmp.chrom_names_and_count
        awk -v OFS='\t' '{print $1+1, $2, $3}' parallel.confirmed > parallel.plusone.confirmed
        awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count parallel.plusone.confirmed > ${SAMPLE}.confirmedsplitreads.bed

        cat lengthfiltered.merged.splitreads.${sample}.bed mapq0_single_unique_choices.bed | sort -k1,1 -k2,2n | uniq > ${sample}.mapq0_and_uniq_splitreads.allfalsepositives.bed
    cd ..
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile