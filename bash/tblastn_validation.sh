while getopts l:g:e:p:q:c:d: option
do
case "${option}"
in
l) LOST_OG=${OPTARG};;
g) LOST_GENOME=${OPTARG};;
e) E_VALUE=${OPTARG};;
p) PIDENT=${OPTARG};;
q) QUERY_COV=${OPTARG};;
c) HIT_COUNT=${OPTARG};;
d) BLAST_DB=${OPTARG};;
esac
done

genome_base=$(basename ${LOST_GENOME})
og_base=$(basename ${LOST_OG})

tblastn -query ${LOST_OG} -subject ${LOST_GENOME} \
    -max_intron_length 3000 \
    -outfmt "6 qacc sacc evalue qlen qstart qend sstart send nident mismatch"  \
    -max_target_seqs 1 > tblastn_${genome_base}_${og_base}

python /global/home/users/pierrj/git/python/parse_tblastn_hits.py tblastn_${genome_base}_${og_base} ${E_VALUE} ${PIDENT} ${QUERY_COV} ${HIT_COUNT} ${genome_base} ${og_base}

if [ -f "${genome_base}_${og_base}.gff3" ]; then
    agat_sp_extract_sequences.pl --gff ${genome_base}_${og_base}.gff3 -f ${LOST_GENOME} -p -o ${genome_base}_${og_base}.fasta &> /dev/null
    blastp -db ${BLAST_DB} -query ${genome_base}_${og_base}.fasta \
        -outfmt "6 qacc sacc evalue qlen qstart qend sstart send nident mismatch"  \
        -max_target_seqs 100 \
        -max_hsps 1 > blastp_${genome_base}_${og_base}
    python /global/home/users/pierrj/git/python/parse_blastp_hits.py blastp_${genome_base}_${og_base} ${og_base} ${genome_base}
fi