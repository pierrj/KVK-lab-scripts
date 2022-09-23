

genome=GCF_000001405.25_GRCh37.p13_genomic.fna
reads_1=empty.fastq
reads_2=empty.fastq

if [ -f "jobqueue" ]; then
    rm jobqueue
fi

while read sample; do 
    echo "python ecc_finder.py map-sr ${genome} ${reads_1} ${reads_2} -r ${genome} -t 1 -o ecc_finder_${sample} > ecc_finder_${sample}/run.out " >> jobqueue
done < mapfile_human_samples

conda activate ecc_finder

nohup parallel -j 32 < jobqueue &