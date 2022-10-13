
cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap

paste genomes_mapfile_no_mgrisae proteomes_mapfile_no_mgrisae > genomes_proteomes_mapfile_no_mgrisae


while read genome proteome; do
    echo $genome
    echo $proteome
    sbatch --export=ALL,genome=$genome,proteome=$proteome /global/home/users/pierrj/git/slurm/pfam_scan_genome.slurm
done < genomes_proteomes_mapfile_no_mgrisae

while read genome proteome; do
    echo $genome
    cp ${genome}_pfam_scan/${genome}.pfamscan.kparse.out all_pfam_scan/${genome}.pfamscan.kparse.out
done < genomes_mapfile_no_mgrisae


cd /global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all

paste wheat_blast_busco_greater_than_90_annotated proteomes_mapfile_no_mgrisae > genomes_proteomes_mapfile_no_mgrisae

while read genome proteome; do
    echo $genome
    echo $proteome
    sbatch --export=ALL,genome=$genome,proteome=$proteome /global/home/users/pierrj/git/slurm/pfam_scan_genome.slurm
done < genomes_proteomes_mapfile_no_mgrisae