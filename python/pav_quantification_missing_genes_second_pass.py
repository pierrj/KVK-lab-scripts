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
from itertools import filterfalse

## USAGE ##
# this script takes an output from orthofinder showing which orthologs in query genomes, the orthologs in a reference genome corresponds to
# if a gene is missing, check it against a list of blastn validated missing genes
# finally output number of genomes the orthogroup is absent from for each gene
# "directory_orthogroups" - directory of orthogroup comparison results
# "directory_missinggenes" - list of fasta/blastn validated missing genes
# "fungap_out" - gene gff file for reference
# "output" - name of output of gene names and number of genomes that are missing the gene

directory_orthogroups = str(sys.argv[1])
directory_missinggenes = str(sys.argv[2])
fungap_out = str(sys.argv[3])
output = str(sys.argv[4])

list_of_files = listdir(directory_orthogroups)
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
    with open(directory_orthogroups + '/' + list_of_files[i], newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        next(file_reader)
        for row in file_reader:
            gene = row[1]
            if len(gene.split(', ')) > 1:
                for ii in range(len(gene.split(', '))): # split cells which lists genes separated by commas
                    split_gene = gene.split(', ')[ii]
                    found_genes_list.append(split_gene)
            else:
                found_genes_list.append(gene)
    genomes_dict[genome] = list(set(guy11_gene_list) - set(found_genes_list))

list_of_files = listdir(directory_missinggenes)

# read in list of genes that are actually missing from the genomes (blastn validated)
genome_dict_actually_missing = {}
for i in range(len(list_of_files)):
    genome = list_of_files[i][:-18]
    putative_missing_list = genomes_dict[genome]
    genome_dict_actually_missing[genome] = []
    not_missing_list = []
    with open(directory_missinggenes + '/' + list_of_files[i], newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        for row in file_reader:
            not_missing_list.append(row[0])
    genome_dict_actually_missing[genome] = list(set(putative_missing_list) - set(not_missing_list))

# count number of genomes each gene is absent from and store as dictionary
pav_final = {}
for i in range(len(guy11_gene_list)):
    gene = guy11_gene_list[i]
    pav_final[gene] = 162
for genome in genome_dict_actually_missing.keys():
    for i in range(len(genome_dict_actually_missing[genome])):
        pav_final[genome_dict_actually_missing[genome][i]] -= 1

# output info stored in genes and pav count dictionary
with open(output, 'w', newline = '') as output_csv:
    w = csv.writer(output_csv, delimiter = '\t')
    for key in pav_final.keys():
        row = [key, pav_final[key]]
        w.writerow(row)