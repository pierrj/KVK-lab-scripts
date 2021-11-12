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
#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
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


module load repeatmasker/4.1.1

cd /global/scratch/users/pierrj/references/te_annotations/rice

# get te annotation file
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_rm.out.gz

gunzip GCF_001433935.1_IRGSP-1.0_rm.out.gz

# filter for copia and gypsy locations
awk -v OFS='\t' '{if ($11 ~ /Copia/) {print $5, $6, $7, $10}}' GCF_001433935.1_IRGSP-1.0_rm.out > rice_copia_locs.bed

awk -v OFS='\t' '{if ($11 ~ /Gypsy/) {print $5, $6, $7, $10}}' GCF_001433935.1_IRGSP-1.0_rm.out > rice_gypsy_locs.bed




cd /global/scratch/users/pierrj/references/te_annotations/arabidopsis

# get te sequences
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_TE.fas

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' TAIR10_TE.fas > TAIR10_TE.singleline.fas

# filter for copia and gypsy names
grep 'Copia' TAIR10_TE.singleline.fas | awk -v FS='|' '{print $1}' > copia_names

grep 'Gypsy' TAIR10_TE.singleline.fas | awk -v FS='|' '{print $1}' > gypsy_names

if [ -f arabidopsis_copia_ltr_seq.fasta ] ; then
    rm arabidopsis_copia_ltr_seq.fasta
fi

# get sequences from names
while read copia;
do
grep -A 1 ${copia} TAIR10_TE.singleline.fas >> arabidopsis_copia_ltr_seq.fasta
done < copia_names

if [ -f arabidopsis_gypsy_ltr_seq.fasta ] ; then
    rm arabidopsis_gypsy_ltr_seq.fasta
fi

# get sequences from names
while read gypsy;
do
grep -A 1 ${gypsy} TAIR10_TE.singleline.fas >> arabidopsis_gypsy_ltr_seq.fasta
done < gypsy_names

RM_LIB=arabidopsis_copia_ltr_seq.fasta
REFERENCE=TAIR10.fasta

# repeatmasking, higher than default cutoff, no small repeats, only copia sequences
RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > arabidopsis_copia_locs.bed

RM_LIB=arabidopsis_gypsy_ltr_seq.fasta
REFERENCE=TAIR10.fasta

# repeatmasking, higher than default cutoff, no small repeats, only gypsy sequences
RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > arabidopsis_gypsy_locs.bed



cd /global/scratch/users/pierrj/references/te_annotations/yeast

# get te annotation file
wget http://sgd-archive.yeastgenome.org/sequence/S288C_reference/other_features/other_features_genomic.fasta.gz

gunzip other_features_genomic.fasta.gz

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' other_features_genomic.fasta > other_features_genomic.singleline.fasta

# gypsy names
grep 'Ty3' other_features_genomic.fasta | awk '{print $1}' > gypsy_names

if [ -f copia_names ] ; then
    rm copia_names
fi

# copia names
grep 'Ty1' other_features_genomic.fasta | awk '{print $1}' >> copia_names
grep 'Ty2' other_features_genomic.fasta | awk '{print $1}' >> copia_names
grep 'Ty4' other_features_genomic.fasta | awk '{print $1}' >> copia_names
grep 'Ty5' other_features_genomic.fasta | awk '{print $1}' >> copia_names

if [ -f yeast_copia_ltr_seq.fasta ] ; then
    rm yeast_copia_ltr_seq.fasta
fi
# get sequences from names
while read copia;
do
grep -A 1 ${copia} other_features_genomic.singleline.fasta >> yeast_copia_ltr_seq.fasta
done < copia_names

if [ -f yeast_gypsy_ltr_seq.fasta ] ; then
    rm yeast_gypsy_ltr_seq.fasta
fi
# get sequences from names
while read gypsy;
do
grep -A 1 ${gypsy} other_features_genomic.singleline.fasta >> yeast_gypsy_ltr_seq.fasta
done < gypsy_names

RM_LIB=yeast_copia_ltr_seq.fasta
REFERENCE=Scer_S288C.fasta
# repeatmasking, higher than default cutoff, no small repeats, only copia sequences
RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > yeast_copia_locs.bed

RM_LIB=yeast_gypsy_ltr_seq.fasta
REFERENCE=Scer_S288C.fasta
# repeatmasking, higher than default cutoff, no small repeats, only gypsy sequences
RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > yeast_gypsy_locs.bed