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
import csv
from os import listdir
import sys

## USAGE ##
# this script takes an output from orthofinder showing which orthologs in query genomes, the orthologs in a reference genome corresponds to
# if a gene is missing, the fasta sequence of that gene is output for blast validation
# "directory" - directory of orthogroup comparison results (i.e. guy11 vs all other genomes)
# "fungap_out" - gene gff file for reference
# "gene_cds" - cds fasta sequences for genes 


directory = str(sys.argv[1])
fungap_out = str(sys.argv[2])
gene_cds = str(sys.argv[3])
list_of_files = listdir(directory)
# remove ncrassa and duplicate guy11 genome
list_of_files.remove('guy11_fungap_out_prot__v__Neurospora_crassa.NC12.pep.all.tsv')
list_of_files.remove('guy11_fungap_out_prot__v__GCA_002368485.1_ASM236848v1_fungap_out_prot.tsv')

# read in list of genes, but only gene locations
guy11_gene_list = []
with open(fungap_out, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        try:
            if row[2] == 'gene':
                guy11_gene_list.append(row[8][3:13]+'.t1')
        except IndexError:
            pass

## read in orthofinder outputs and organize them into a dictionary with teh names of the genomes as the keys
genomes_dict = {}
for i in range(len(list_of_files)):
    if list_of_files[i] == 'guy11_fungap_out_prot__v__Neurospora_crassa.NC12.pep.all.tsv':
        continue
    genome = list_of_files[i][26:-20]
    found_genes_list = []
    with open(directory + '/' + list_of_files[i], newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        next(file_reader) # skip first line
        for row in file_reader:
            gene = row[1]
            if len(gene.split(', ')) > 1:
                for ii in range(len(gene.split(', '))): # split cells which lists genes separated by commas
                    split_gene = gene.split(', ')[ii]
                    found_genes_list.append(split_gene)
            else:
                found_genes_list.append(gene)
    genomes_dict[genome] = list(set(guy11_gene_list) - set(found_genes_list))

# open gene cds and store as dictionary
gene_fastas = {}
with open(gene_cds, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        gene = row[0][1:]
        row2 = next(file_reader)
        fasta = row2[0]
        gene_fastas[gene] = fasta

# output fasta file 
for key in genomes_dict.keys():
    with open(key+'_missing_genes.fasta', 'w', newline = '') as output_csv:
        w = csv.writer(output_csv, delimiter = '\t')
        for i in range(len(genomes_dict[key])):
            gene = genomes_dict[key][i]
            gene_list = ['>' + genomes_dict[key][i]]
            fasta = [gene_fastas[gene]]
            w.writerow(gene_list)
            w.writerow(fasta)