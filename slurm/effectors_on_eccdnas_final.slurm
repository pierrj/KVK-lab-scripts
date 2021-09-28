#!/bin/bash
#SBATCH --job-name=effectors_on_eccdnas
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


## predict candidate effectors

# bin/signalp -fasta /Users/pierrj/references/annotated_rice_genomes_12_29_20/guy11_fungap_out_prot.faa -prefix signalp_guy11

# awk '{if ($2 == "SP(Sec/SPI)") {print $1}}' signalp_guy11_summary.signalp5 > signalp_secrete_points_names

# awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' /Users/pierrj/references/annotated_rice_genomes_12_29_20/guy11_fungap_out_prot.faa > guy11_fungap_out_prot.singleline.faa
# while read gene; do
# grep -A1 ${gene} guy11_fungap_out_prot.singleline.faa >> guy11_signalp_proteins.faa
# done < signalp_secrete_points_names

# ./tmhmm guy11_signalp_proteins.faa > tmhmm_output

# grep "Number" tmhmm_output | awk '{if ($7 == 0){print $2}} ' > guy11_signalp_notmhmm_protein_names

# while read gene; do
# grep -A1 ${gene} guy11_fungap_out_prot.singleline.faa >> guy11_signalp_notmhmm_proteins.faa
# done < guy11_signalp_notmhmm_protein_names

# python EffectorP.py -i guy11_signalp_notmhmm_proteins.faa

# grep "Effector probability" effectorp_output | awk '{print $1}' > guy11_effector_protein_names

GENOME_FASTA=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
SAMPLE=G3
SRA_LIST=/global/scratch/pierrj/eccDNA/pipeline_tests/eccdna_gene_expression/zhang_et_al_data/zhang_et_al_rnaseq_accessions.txt
THREADS=${SLURM_NTASKS}
LIB_TYPE=1
GFF_FILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.gff3
SAMPLE_MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/rnaseq_comparison_wrapper/anyoverlap_normalization/G3.sample_mapfile
GENE_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.renamed.bed
OUTPUT_NAME=G3

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common

ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/neverfound/G3.all_eccdnas_mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do
cd ${sample}
    realpath ${sample}.ecc_caller_out.splitreads.bed >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/common

/global/home/users/pierrj/git/bash/output_neverfound_common_genes.sh -g ${GENOME_FASTA} \
    -s ${SAMPLE} -a ${SRA_LIST} -t ${THREADS} -l ${LIB_TYPE} -f ${GFF_FILE} \
    -e ${ECCDNA_MAPFILE} -m ${SAMPLE_MAPFILE} -n n

## analyze the split reads for effectors and noneffectors