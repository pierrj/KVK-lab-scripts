import csv
import pandas as pd
import numpy as np
import sys

orthogroup_tsv = sys.argv[1]
guy11_data_dir = sys.argv[2]

og_w_guy11_df = pd.read_csv(orthogroup_tsv, dtype='string', sep='\t', index_col = 0)

## to get the og that a gene belongs to for orthogrouping with guy11
df_lol = og_w_guy11_df.values.tolist()

og_dict_w_guy11 = {}

for i, row in enumerate(df_lol):
    og = og_w_guy11_df.index[i]
    for cell in row:
        if not pd.isnull(cell):
            for protein in cell.split(', '):
                og_dict_w_guy11[protein] = og

## to get all genes associated with one OG for orthogrouping with guy11
genes_per_og_w_guy11 = {}

for gene in og_dict_w_guy11:
    og = og_dict_w_guy11[gene]
    if og not in genes_per_og_w_guy11:
        genes_per_og_w_guy11[og] = []
    genes_per_og_w_guy11[og].append(gene)


input_files = [
    'guy11_H3K27ac_per_gene.txt',
    'guy11_H3K27me3_per_gene.txt',
    'guy11_H3K36me3_per_gene.txt',
    'guy11_zhang_et_al_2019_complete_medium_expression.txt',
    'guy11_zhang_et_al_2019_in_planta_expression.txt',
    'guy11_eccdnacov_per_gene.txt',
    'guy11_methylation_per_gene.txt'
]

output_files =[
    'H3K27ac_per_og.txt',
    'H3K27me3_per_og.txt',
    'H3K36me3_per_og.txt',
    'zhang_et_al_2019_complete_medium_expression_per_og.txt',
    'zhang_et_al_2019_in_planta_expression_per_og.txt',
    'eccdnacov_per_og.txt',
    'methylation_per_og.txt'
]


for i in range(len(input_files)):
    input_file = input_files[i]
    print(input_file)
    output_file = output_files[i]
    signal_per_gene_dict = {}
    with open(guy11_data_dir+'/'+input_file, newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        for row in file_reader:
            if "ID=" in row[0]:
                gene = row[0].split(';')[0][3:]
            else:
                gene = row[0][:-2]
            signal_per_gene_dict[gene] = float(row[1])
    for og in genes_per_og_w_guy11: ##this is specifically for genes that have no cytosines that can be methylated
        for gene in genes_per_og_w_guy11[og]:
            if "GUY11" in gene and gene not in signal_per_gene_dict:
                signal_per_gene_dict[gene] = 0
    og_signal_w_guy11 = {}
    for og in genes_per_og_w_guy11:
        og_signal_w_guy11[og] = []
        for gene in genes_per_og_w_guy11[og]:
            if "GUY11" in gene:
                og_signal_w_guy11[og].append(signal_per_gene_dict[gene])
    og_signal_w_guy11_averaged = {}
    for og in og_signal_w_guy11:
        lst = og_signal_w_guy11[og]
        if len(lst) > 1:
            og_signal_w_guy11_averaged[og] = sum(lst) / len(lst)
        elif len(lst) == 1:
            og_signal_w_guy11_averaged[og] = lst[0]
        elif len(lst) == 0:
            pass
        else:
            print('wtf')
    median_value = np.median(list(og_signal_w_guy11_averaged.values()))
    imputed_values = []
    for og in og_signal_w_guy11:
        if og not in og_signal_w_guy11_averaged:
            imputed_values.append(og)
            og_signal_w_guy11_averaged[og] = median_value
    print(len(imputed_values))
    with open(output_file, 'w', newline = '') as output_csv:
        w = csv.writer(output_csv, delimiter = '\t')
        for key in og_signal_w_guy11_averaged:
            w.writerow([key, og_signal_w_guy11_averaged[key]])