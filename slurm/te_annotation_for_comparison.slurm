#!/bin/bash
#SBATCH --job-name=te_annotation_for_comparison
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

module load repeatmasker/4.1.1

cd /global/scratch/users/pierrj/references/te_annotations/rice

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_rm.out.gz

gunzip GCF_001433935.1_IRGSP-1.0_rm.out.gz

awk -v OFS='\t' '{if ($11 ~ /Copia/) {print $5, $6, $7, $10}}' GCF_001433935.1_IRGSP-1.0_rm.out > rice_copia_locs.bed

awk -v OFS='\t' '{if ($11 ~ /Gypsy/) {print $5, $6, $7, $10}}' GCF_001433935.1_IRGSP-1.0_rm.out > rice_gypsy_locs.bed




cd /global/scratch/users/pierrj/references/te_annotations/arabidopsis

wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_TE.fas

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' TAIR10_TE.fas > TAIR10_TE.singleline.fas

grep 'Copia' TAIR10_TE.singleline.fas | awk -v FS='|' '{print $1}' > copia_names

grep 'Gypsy' TAIR10_TE.singleline.fas | awk -v FS='|' '{print $1}' > gypsy_names

if [ -f arabidopsis_copia_ltr_seq.fasta ] ; then
    rm arabidopsis_copia_ltr_seq.fasta
fi

while read copia;
do
grep -A 1 ${copia} TAIR10_TE.singleline.fas >> arabidopsis_copia_ltr_seq.fasta
done < copia_names

if [ -f arabidopsis_gypsy_ltr_seq.fasta ] ; then
    rm arabidopsis_gypsy_ltr_seq.fasta
fi

while read gypsy;
do
grep -A 1 ${gypsy} TAIR10_TE.singleline.fas >> arabidopsis_gypsy_ltr_seq.fasta
done < gypsy_names

RM_LIB=arabidopsis_copia_ltr_seq.fasta
REFERENCE=TAIR10.fasta

RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > arabidopsis_copia_locs.bed

RM_LIB=arabidopsis_gypsy_ltr_seq.fasta
REFERENCE=TAIR10.fasta

RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > arabidopsis_gypsy_locs.bed



cd /global/scratch/users/pierrj/references/te_annotations/yeast

wget http://sgd-archive.yeastgenome.org/sequence/S288C_reference/other_features/other_features_genomic.fasta.gz

gunzip other_features_genomic.fasta.gz

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' other_features_genomic.fasta > other_features_genomic.singleline.fasta

grep 'Ty3' other_features_genomic.fasta | awk '{print $1}' > gypsy_names

if [ -f copia_names ] ; then
    rm copia_names
fi

grep 'Ty1' other_features_genomic.fasta | awk '{print $1}' >> copia_names
grep 'Ty2' other_features_genomic.fasta | awk '{print $1}' >> copia_names
grep 'Ty4' other_features_genomic.fasta | awk '{print $1}' >> copia_names
grep 'Ty5' other_features_genomic.fasta | awk '{print $1}' >> copia_names

if [ -f yeast_copia_ltr_seq.fasta ] ; then
    rm yeast_copia_ltr_seq.fasta
fi

while read copia;
do
grep -A 1 ${copia} other_features_genomic.singleline.fasta >> yeast_copia_ltr_seq.fasta
done < copia_names

if [ -f yeast_gypsy_ltr_seq.fasta ] ; then
    rm yeast_gypsy_ltr_seq.fasta
fi

while read gypsy;
do
grep -A 1 ${gypsy} other_features_genomic.singleline.fasta >> yeast_gypsy_ltr_seq.fasta
done < gypsy_names

RM_LIB=yeast_copia_ltr_seq.fasta
REFERENCE=Scer_S288C.fasta

RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > yeast_copia_locs.bed

RM_LIB=yeast_gypsy_ltr_seq.fasta
REFERENCE=Scer_S288C.fasta

RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > yeast_gypsy_locs.bed


## for moryzae, probably use wickersoft called LTR-TEs I think

cd /global/scratch/users/pierrj/references/te_annotations/moryzae

cp /global/scratch/pierrj/references/moryzae_ltr_elements.fasta .

grep 'copia' moryzae_ltr_elements.fasta > copia_names

grep -v 'copia' moryzae_ltr_elements.fasta | grep '>' > gypsy_names

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' moryzae_ltr_elements.fasta > moryzae_ltr_elements.singleline.fasta

if [ -f moryzae_copia_ltr_seq.fasta ] ; then
    rm moryzae_copia_ltr_seq.fasta
fi

while read copia;
do
grep -A 1 ${copia} moryzae_ltr_elements.singleline.fasta >> moryzae_copia_ltr_seq.fasta
done < copia_names

if [ -f moryzae_gypsy_ltr_seq.fasta ] ; then
    rm moryzae_gypsy_ltr_seq.fasta
fi

while read gypsy;
do
grep -A 1 ${gypsy} moryzae_ltr_elements.singleline.fasta >> moryzae_gypsy_ltr_seq.fasta
done < gypsy_names

RM_LIB=moryzae_copia_ltr_seq.fasta
REFERENCE=guy11_genome_baoetal2017.fasta

RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > moryzae_copia_locs.bed

RM_LIB=moryzae_gypsy_ltr_seq.fasta
REFERENCE=guy11_genome_baoetal2017.fasta

RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > moryzae_gypsy_locs.bed