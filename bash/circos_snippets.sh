
NEVERFOUND_FILE=/global/scratch/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/G3.neverfound.genes
COMMON_FILE=/global/scratch/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/G3.common.genes

if [ -f "neverfound_locs" ]; then
    rm neverfound_locs
fi

while read gene; do
    gene_trimmed=$(echo $gene | awk '{print substr($1, 0, length($1)-3)}')
    grep ${gene_trimmed} ${GENE_BEDFILE} >> neverfound_locs
done < ${NEVERFOUND_FILE}

if [ -f "common_locs" ]; then
    rm common_locs
fi

while read gene; do
    gene_trimmed=$(echo $gene | awk '{print substr($1, 0, length($1)-3)}')
    grep ${gene_trimmed} ${GENE_BEDFILE} >> common_locs
done < ${COMMON_FILE}

cat neverfound_locs common_locs ${GENE_BEDFILE} | sort | uniq -u > other_locs


awk -v OFS='\t' '{print $1, $2, $3, 0}' neverfound_locs > neverfound_locs_color
awk -v OFS='\t' '{print $1, $2, $3, 1}' other_locs > other_locs_color
awk -v OFS='\t' '{print $1, $2, $3, 2}' common_locs > common_locs_color

cat neverfound_locs_color other_locs_color common_locs_color | sort -k1,1 -k2,2n | \
    awk -v OFS='\t' '{print substr($1,11,2), $2, $3, $4}' | | awk '{gsub ("^0*", "", $0); gsub ("/0*", "/", $0); print}' > gene_status