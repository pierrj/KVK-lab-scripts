import csv
from os import listdir
import sys
import collections

directory = str(sys.argv[1])
outgroup_file = str(sys.argv[2])
output = str(sys.argv[3])

list_of_files = os.listdir(directory)

genes_dict = {}
for i in range(len(list_of_files)):
    print(list_of_files[i])
    if list_of_files[i] == outgroup_file: ## skip outgroup
        continue
    with open(directory + '/' + list_of_files[i], newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        next(file_reader)
        for row in file_reader:
            gene = row[1]
            if len(key_trimmed.split(', ')) > 1:
                for ii in range(len(gene.split(', '))):
                    split_gene = gene.split(', ')[ii]
                    if split_gene not in genes_dict:
                        genes_dict[split_gene] = 1
                    else:
                        genes_dict[split_gene] += 1
            else:
                if gene not in genes_dict:
                    genes_dict[gene] = 1
                else:
                    genes_dict[gene] += 1

genes_list = []
for key, val in genes_dict.items():
    genes_list.append([key, val])

genes_list_sorted = sorted(genes_list, key=lambda x: x[0])

output = 'pav_per_gene' # str(sys.argv[2])
with open(output, 'w', newline = '') as output_csv:
    w = csv.writer(output_csv, delimiter = '\t')
    for i in range(len(genes_list_sorted)):
        row = [genes_list_sorted[i][0], genes_list_sorted[i][1]]
        w.writerow(row)