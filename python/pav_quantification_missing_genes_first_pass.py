import csv
from os import listdir
import sys


directory = str(sys.argv[1])
fungap_out = str(sys.argv[2])
gene_cds = str(sys.argv[3])
list_of_files = listdir(directory)
list_of_files.remove('guy11_fungap_out_prot__v__Neurospora_crassa.NC12.pep.all.tsv')
list_of_files.remove('guy11_fungap_out_prot__v__GCA_002368485.1_ASM236848v1_fungap_out_prot.tsv')

guy11_gene_list = []
with open(fungap_out, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        try:
            if row[2] == 'gene':
                guy11_gene_list.append(row[8][3:13]+'.t1')
        except IndexError:
            pass

genomes_dict = {}
for i in range(len(list_of_files)):
    if list_of_files[i] == 'guy11_fungap_out_prot__v__Neurospora_crassa.NC12.pep.all.tsv':
        continue
    genome = list_of_files[i][26:-20]
    found_genes_list = []
    with open(directory + '/' + list_of_files[i], newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        next(file_reader)
        for row in file_reader:
            gene = row[1]
            if len(gene.split(', ')) > 1:
                for ii in range(len(gene.split(', '))):
                    split_gene = gene.split(', ')[ii]
                    found_genes_list.append(split_gene)
            else:
                found_genes_list.append(gene)
    genomes_dict[genome] = list(set(guy11_gene_list) - set(found_genes_list))

gene_fastas = {}
with open(gene_cds, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        gene = row[0][1:]
        row2 = next(file_reader)
        fasta = row2[0]
        gene_fastas[gene] = fasta

for key in genomes_dict.keys():
    with open(key+'_missing_genes.fasta', 'w', newline = '') as output_csv:
        w = csv.writer(output_csv, delimiter = '\t')
        for i in range(len(genomes_dict[key])):
            gene = genomes_dict[key][i]
            gene_list = ['>' + genomes_dict[key][i]]
            fasta = [gene_fastas[gene]]
            w.writerow(gene_list)
            w.writerow(fasta)