


for i in {1..1000}; do 
    bedtools shuffle -i noduplicates.all.nonltr_eccs.txt -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes -excl /global/scratch/pierrj/references/te_annotations/moryzae/moryzae_copiaandgypsy_locs.bed > shuffled
    python /global/home/users/pierrj/git/python/check_overlap_ecc_and_regions.py zhang_all.ciri_out shuffled 100
done


for i in {1..1000}; do 
    bedtools shuffle -i all.nonltr_splitreads.txt -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes -excl /global/scratch/pierrj/references/te_annotations/moryzae/moryzae_copiaandgypsy_locs.bed > shuffled
    python /global/home/users/pierrj/git/python/check_overlap_ecc_and_regions.py zhang_all.ciri_out shuffled 100
done


for i in {1..1000}; do 
    bedtools shuffle -i all.nonltr_splitreads.txt -g /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes -excl /global/scratch/pierrj/references/te_annotations/moryzae/moryzae_copiaandgypsy_locs.bed > shuffled
    python /global/home/users/pierrj/git/python/check_hotspot_ecc_and_regions.py zhang_all.ciri_out shuffled
done