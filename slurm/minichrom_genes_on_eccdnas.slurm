#!/bin/bash
#SBATCH --job-name=minichrom_genes_on_eccdnas
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/

module load STAR

cp /global/scratch/pierrj/eccDNA/pipeline_tests/minichromosomes/GCA_900474655.3_FR13_genome_polished_no_mito_genomic.fna .
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' GCA_900474655.3_FR13_genome_polished_no_mito_genomic.fna > GCA_900474655.3_FR13_genome_polished_no_mito_genomic.fna.singleline
echo -e UEMA03000009.1'\n'UEMA03000023.1 > FR13_minichrom_scaffolds
if [ -f FR13_minichromosome.fasta ] ; then
    rm FR13_minichromosome.fasta
fi

while read scaffold;
do
grep -A 1 ${scaffold} GCA_900474655.3_FR13_genome_polished_no_mito_genomic.fna.singleline >> FR13_minichromosome.fasta
done < FR13_minichrom_scaffolds

cp /global/scratch/pierrj/eccDNA/pipeline_tests/minichromosomes/GCA_900474475.3_CD156_genome_polished_no_mito_genomic.fna .
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' GCA_900474475.3_CD156_genome_polished_no_mito_genomic.fna > GCA_900474475.3_CD156_genome_polished_no_mito_genomic.fna.singleline
echo -e UELZ03000011.1"\n"UELZ03000012.1"\n"UELZ03000014.1 > CD156_minichrom_scaffolds

if [ -f CD156_minichromosome.fasta ] ; then
    rm CD156_minichromosome.fasta
fi

while read scaffold;
do
grep -A 1 ${scaffold} GCA_900474475.3_CD156_genome_polished_no_mito_genomic.fna.singleline >> CD156_minichromosome.fasta
done < CD156_minichrom_scaffolds

cp /global/scratch/pierrj/eccDNA/pipeline_tests/minichromosomes/GCA_900474175.3_US71_genome_polished_no_mito_genomic.fna .
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' GCA_900474175.3_US71_genome_polished_no_mito_genomic.fna > GCA_900474175.3_US71_genome_polished_no_mito_genomic.fna.singleline
echo -e UCNY03000014.1"\n"UCNY03000021.1"\n"UCNY03000028.1"\n"UCNY03000029.1 > US71_minichrom_scaffolds

if [ -f US71_minichromosome.fasta ] ; then
    rm US71_minichromosome.fasta
fi

while read scaffold;
do
grep -A 1 ${scaffold} GCA_900474175.3_US71_genome_polished_no_mito_genomic.fna.singleline >> US71_minichromosome.fasta
done < US71_minichrom_scaffolds


CDS_FASTA=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20_prot.faa

module load exonerate

isolate=FR13

exonerate --model protein2genome ${CDS_FASTA} ${isolate}_minichromosome.fasta --showtargetgff yes --bestn 1 --showalignment no --verbose 0 --showvulgar no > ${isolate}_genes.gff

awk '$3=="gene"' ${isolate}_genes.gff | awk '$19>70' | awk '{print $13}' > ${isolate}_genes_on_minichromosome

isolate=CD156

exonerate --model protein2genome ${CDS_FASTA} ${isolate}_minichromosome.fasta --showtargetgff yes --bestn 1 --showalignment no --verbose 0 --showvulgar no > ${isolate}_genes.gff

awk '$3=="gene"' ${isolate}_genes.gff | awk '$19>70' | awk '{print $13}' > ${isolate}_genes_on_minichromosome

isolate=US71

exonerate --model protein2genome ${CDS_FASTA} ${isolate}_minichromosome.fasta --showtargetgff yes --bestn 1 --showalignment no --verbose 0 --showvulgar no > ${isolate}_genes.gff

awk '$3=="gene"' ${isolate}_genes.gff | awk '$19>70' | awk '{print $13}' > ${isolate}_genes_on_minichromosome


cat FR13_genes_on_minichromosome CD156_genes_on_minichromosome US71_genes_on_minichromosome > all_genes_on_minichromosome


ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/G3.all_eccdnas_mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do
cd ${sample}
    realpath ${sample}.ecc_caller_out.splitreads.bed >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile

GENOME_FASTA=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
SAMPLE=G3
SRA_LIST=/global/scratch/pierrj/eccDNA/pipeline_tests/eccdna_gene_expression/zhang_et_al_data/zhang_et_al_rnaseq_accessions.txt
THREADS=${SLURM_NTASKS}
LIB_TYPE=1
GFF_FILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.gff3
SAMPLE_MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/rnaseq_comparison_wrapper/anyoverlap_normalization/G3.sample_mapfile
GENE_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.renamed.bed
OUTPUT_NAME=G3

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common/

/global/home/users/pierrj/git/bash/output_neverfound_common_genes.sh -g ${GENOME_FASTA} \
    -s ${SAMPLE} -a ${SRA_LIST} -t ${THREADS} -l ${LIB_TYPE} -f ${GFF_FILE} \
    -e ${ECCDNA_MAPFILE} -m ${SAMPLE_MAPFILE} -n n

## just use the splitread counts again for analysis here