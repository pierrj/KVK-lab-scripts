while getopts l:g:e:p:q:c:f:b: option
do
case "${option}"
in
l) LOST_OG=${OPTARG};;
g) LOST_GENOME=${OPTARG};;
e) E_VALUE=${OPTARG};;
p) PIDENT=${OPTARG};;
q) QUERY_COV=${OPTARG};;
c) HIT_COUNT=${OPTARG};;
f) GENE_GFF=${OPTARG};;
b) GENE_OVERLAP=${OPTARG};;
esac
done

genome_base=$(basename ${LOST_GENOME})
og_base=$(basename ${LOST_OG})

tblastn -query ${LOST_OG} -subject ${LOST_GENOME} \
    -max_intron_length 3000 \
    -outfmt "6 qacc sacc evalue qlen qstart qend sstart send nident mismatch"  \
    -max_target_seqs 1 > tblastn_${genome_base}_${og_base}

if [ ! -f "${GENE_GFF}.justgenes" ]; then
    awk '$3 == "gene"' ${GENE_GFF} > ${GENE_GFF}.justgenes
fi

python /global/home/users/pierrj/git/python/parse_tblastn_hits.py tblastn_${genome_base}_${og_base} ${E_VALUE} ${PIDENT} ${QUERY_COV} ${HIT_COUNT} ${genome_base} ${og_base} ${GENE_GFF}.justgenes ${GENE_OVERLAP}