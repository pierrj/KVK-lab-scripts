#!/bin/bash
#SBATCH --job-name=guy11_to_70-15_names
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

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/gene_names

query=guy11_fungap_out_12_28_20_prot.faa

# blast annotated guy11 proteome to 70-15 uniprot proteome
blastp -query ${query} -db uniport_moryzae_70-15.fasta \
    -out ${query}.out -num_threads 20 -outfmt "6 qacc sacc evalue qstart qend sstart send qlen slen pident length stitle " -max_target_seqs 1


#80% sequence identity over 80% reciprocal sequence length
# for reciprocal, take length of alignment and make sure it is greater than 80% of both
awk '$10 >= 80 && $11/$8 >= 0.8 && $11/$9 >= 0.8' ${query}.out > ${query}.out.filtered

awk -v OFS='\t' '{print $1, substr($(NF-2),4)}' ${query}.out.filtered > ${query}.out.filtered.names

awk -v OFS='\t' '{print $1}' ${query}.out.filtered > annotated_genes

grep '>' $query | awk '{print substr($1, 2)}' > all_genes

# add unknowns to those that weren't identified
cat annotated_genes all_genes | awk -v OFS='\t' '{ print $1, "unknown"}' | sort | uniq -u | cat - ${query}.out.filtered.names | sort -k1,1 | awk -v OFS='\t' '{print $2, $1}' > guy11_70_15_gene_names.txt