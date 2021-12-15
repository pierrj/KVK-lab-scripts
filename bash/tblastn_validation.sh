while getopts l:g:e:p:q:c: option
do
case "${option}"
in
l) LOST_OG=${OPTARG};;
g) LOST_GENOME=${OPTARG};;
e) E_VALUE=${OPTARG};;
p) PIDENT=${OPTARG};;
q) QUERY_COV=${OPTARG};;
c) HIT_COUNT=${OPTARG};;
esac
done

tblastn -query ${LOST_OG} -subject ${LOST_GENOME} \
    -max_intron_length 3000 \
    -outfmt "6 qacc sacc evalue qlen qstart qend sstart send  pident qcovs qcovhsp"  \
    -max_target_seqs 1 > tblastn_${LOST_GENOME}_${LOST_OG}

python parse_tblastn_hits.py tblastn_${LOST_GENOME}_${LOST_OG} ${E_VALUE} ${PIDENT} ${QUERY_COV} ${HIT_COUNT} ${LOST_OG}