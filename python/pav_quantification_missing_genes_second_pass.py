import csv
from os import listdir
import sys
from itertools import filterfalse

directory_orthogroups = str(sys.argv[1])
directory_missinggenes = str(sys.argv[2])
fungap_out = str(sys.argv[3])
output = str(sys.argv[4])

list_of_files = listdir(directory_orthogroups)
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
    with open(directory_orthogroups + '/' + list_of_files[i], newline = '') as file:
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

list_of_files = listdir(directory_missinggenes)

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

pav_final = {}
for i in range(len(guy11_gene_list)):
    gene = guy11_gene_list[i]
    pav_final[gene] = 162
for genome in genome_dict_actually_missing.keys():
    for i in range(len(genome_dict_actually_missing[genome])):
        pav_final[genome_dict_actually_missing[genome][i]] -= 1

with open(output, 'w', newline = '') as output_csv:
    w = csv.writer(output_csv, delimiter = '\t') ## lineterminator="\n" might fix issues with windows
    for key in pav_final.keys():
        row = [key, pav_final[key]]
        w.writerow(row)