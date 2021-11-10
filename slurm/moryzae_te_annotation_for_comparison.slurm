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
#!/bin/bash
#SBATCH --job-name=moryzae_te_annotation_for_comparison
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/references/te_annotations/moryzae

module load repeatmasker/4.1.1

# from ltr consensus sequences generated using wickersoft
cp /global/scratch/pierrj/references/moryzae_ltr_elements.fasta .

# get names of copia and gypsy sequences
grep 'copia' moryzae_ltr_elements.fasta > copia_names

grep -v 'copia' moryzae_ltr_elements.fasta | grep '>' > gypsy_names

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' moryzae_ltr_elements.fasta > moryzae_ltr_elements.singleline.fasta

if [ -f moryzae_copia_ltr_seq.fasta ] ; then
    rm moryzae_copia_ltr_seq.fasta
fi

# just get sequences associated with names
while read copia;
do
grep -A 1 ${copia} moryzae_ltr_elements.singleline.fasta >> moryzae_copia_ltr_seq.fasta
done < copia_names

if [ -f moryzae_gypsy_ltr_seq.fasta ] ; then
    rm moryzae_gypsy_ltr_seq.fasta
fi

# just get sequences associated with names
while read gypsy;
do
grep -A 1 ${gypsy} moryzae_ltr_elements.singleline.fasta >> moryzae_gypsy_ltr_seq.fasta
done < gypsy_names

RM_LIB=moryzae_copia_ltr_seq.fasta
REFERENCE=guy11_genome_baoetal2017.fasta

# repeatmasking, higher than default cutoff, no small repeats, only copia sequences
RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > moryzae_copia_locs.bed

RM_LIB=moryzae_gypsy_ltr_seq.fasta
REFERENCE=guy11_genome_baoetal2017.fasta

# repeatmasking, higher than default cutoff, no small repeats, only gypsy sequences
RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > moryzae_gypsy_locs.bed

# both gypsy and copia
cat moryzae_gypsy_locs.bed moryzae_copia_locs.bed > moryzae_copiaandgypsy_locs.bed

# from te sequences generated from repeatmodeler
cp /global/scratch/pierrj/eccDNA/pipeline_tests/ltr_ecc_filter/repbase_plus_guy11_modeler_rm_database.fasta .

if [ -f LTR_fasta_names.txt ] ; then
    rm LTR_fasta_names.txt
fi
# awkward filtering for ltr te sequences
awk '{if ($1 ~ /family/ && $1 ~ /LTR/) {print $0}}' repbase_plus_guy11_modeler_rm_database.fasta >> LTR_fasta_names.txt
awk '$2 == "Copia"' repbase_plus_guy11_modeler_rm_database.fasta >> LTR_fasta_names.txt
awk '$2 == "Gypsy"' repbase_plus_guy11_modeler_rm_database.fasta >> LTR_fasta_names.txt
awk '$2=="LTR"' repbase_plus_guy11_modeler_rm_database.fasta >> LTR_fasta_names.txt
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < repbase_plus_guy11_modeler_rm_database.fasta > repbase_plus_guy11_modeler_rm_database.singleline.fasta
## just need first field
mv LTR_fasta_names.txt LTR_fasta_names.txt.old
awk '{print $1}' LTR_fasta_names.txt.old > LTR_fasta_names.txt

grep '>' repbase_plus_guy11_modeler_rm_database.singleline.fasta | awk '{print $1}' > all_fasta_names.txt
comm -23 <(sort all_fasta_names.txt) <(sort LTR_fasta_names.txt) > nonLTR_fasta_names.txt

if [ -f repbase_plus_guy11_modeler_rm_database.nonLTR_only.fasta ] ; then
    rm repbase_plus_guy11_modeler_rm_database.nonLTR_only.fasta
fi
# nonltr tes
while read nonltr; do
grep -A1 "${nonltr}[[:space:]]" repbase_plus_guy11_modeler_rm_database.singleline.fasta >> repbase_plus_guy11_modeler_rm_database.nonLTR_only.fasta
done < nonLTR_fasta_names.txt

module load repeatmasker/4.1.1

RM_LIB=repbase_plus_guy11_modeler_rm_database.nonLTR_only.fasta
REFERENCE=guy11_genome_baoetal2017.fasta

# repeatmasking, higher than default cutoff, no small repeats, only nonltr sequences
RepeatMasker -cutoff 250 -gff -no_is -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > moryzae_non_ltr_repeats.bed

cp /global/scratch/users/pierrj/references/moryzae_ltr_consensus_sequences.fasta .

RM_LIB=moryzae_ltr_consensus_sequences.fasta
REFERENCE=guy11_genome_baoetal2017.fasta

# repeatmasking, higher than default cutoff, no small repeats, only ltr sequences
RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{if ($5-$4>3000 && $6 < 10) {print $1, $4, $5, $10}}' ${REFERENCE}.out.gff > moryzae_high_conf_ltr_tes_for_plotprofile.bed

# from ltr consensus sequences generated using wickersoft again
cp /global/scratch/users/pierrj/references/moryzae_ltr_consensus_sequences.fasta .

RM_LIB=moryzae_ltr_consensus_sequences.fasta
REFERENCE=guy11_genome_baoetal2017.fasta

# repeatmasking, higher than default cutoff, no small repeats, only ltr sequences
RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > moryzae_high_conf_ltr_tes_nolengthfilter.bed

cp /global/scratch/users/pierrj/references/moryzae_ltr_elements.fasta .

awk '{if ($0 ~ /Moryzae_copia1_LTR/ ) {print ">Moryzae_copia1_LTR"} else {print $0}}' moryzae_ltr_elements.fasta | \
    awk '{if ($0 ~ /Moryzae_copia1_I/ ) {print ">Moryzae_copia1_INTERNAL"} else {print $0}}' > moryzae_ltr_tes_split_by_element.fasta

RM_LIB=moryzae_ltr_tes_split_by_element.fasta
REFERENCE=guy11_genome_baoetal2017.fasta

# repeatmasking, higher than default cutoff, no small repeats, only ltr sequences
RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > moryzae_high_conf_ltr_tes_nolengthfilter_individual_elements.bed

# from repeatmodeler
RM_LIB=repbase_plus_guy11_modeler_rm_database.fasta
REFERENCE=guy11_genome_baoetal2017.fasta

# repeatmasking, higher than default cutoff, no small repeats, all te sequences, no simple repeats
RepeatMasker -cutoff 250 -gff -nolow -no_is -norna -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > moryzae_all_tes_norepeats.bed

RM_LIB=repbase_plus_guy11_modeler_rm_database.fasta
REFERENCE=guy11_genome_baoetal2017.fasta

# repeatmasking, higher than default cutoff, no small repeats, all te sequences, with simple repeats
RepeatMasker -cutoff 250 -gff -no_is -pa ${SLURM_NTASKS} -lib ${RM_LIB} -dir . ${REFERENCE}

awk -v OFS='\t' '{print $1, $4, $5, $10}' ${REFERENCE}.out.gff > moryzae_all_tes_wrepeats.bed