#MIT License
#
#Copyright (c) 2023 Pierre Michel Joubert
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

## tblastn lost orthogroup file to genome
tblastn -query ${LOST_OG} -subject ${LOST_GENOME} \
    -max_intron_length 3000 \
    -outfmt "6 qacc sacc evalue qlen qstart qend sstart send nident mismatch"  \
    -max_target_seqs 1 > tblastn_${genome_base}_${og_base}

# parse and filter tblastn hits 
python /global/home/users/pierrj/git/python/parse_tblastn_hits.py tblastn_${genome_base}_${og_base} ${E_VALUE} ${PIDENT} ${QUERY_COV} ${HIT_COUNT} ${genome_base} ${og_base}

# if any hits are left after filtration
if [ -f "${genome_base}_${og_base}.gff3" ]; then
    # get fasta from gff
    agat_sp_extract_sequences.pl --gff ${genome_base}_${og_base}.gff3 -f ${LOST_GENOME} -p -o ${genome_base}_${og_base}.fasta &> /dev/null
    # blast fasta against all orthogroups blast db
    blastp -db ${BLAST_DB} -query ${genome_base}_${og_base}.fasta \
        -outfmt "6 qacc sacc evalue qlen qstart qend sstart send nident mismatch"  \
        -max_target_seqs 100 \
        -max_hsps 1 > blastp_${genome_base}_${og_base}
    # parse and filter final blast p hits
    python /global/home/users/pierrj/git/python/parse_blastp_hits.py blastp_${genome_base}_${og_base} ${og_base} ${genome_base}
fi