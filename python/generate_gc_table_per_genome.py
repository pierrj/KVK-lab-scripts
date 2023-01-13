#MIT License
#
#Copyright (c) 2023 Pierre Michel Joubert
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
import sys

input_gene_gc = sys.argv[1]
input_flanking_gc = sys.argv[2]
genome_name = sys.argv[3]

## read in gene gc content
with open(input_gene_gc, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    # make dictionary
    gene_dict = {}
    for row in file_reader:
        # skip first row
        if "usercol" in row[0]:
            continue
        gene_name = row[3]
        # calc gc content from base counts
        gene_gc_count = int(row[7]) + int(row[8])
        gene_at_count = int(row[6]) + int(row[9])
        gene_gc_perc = gene_gc_count/(gene_gc_count + gene_at_count)
        # add to dict
        gene_dict[gene_name] = []
        gene_dict[gene_name].append(gene_gc_perc)

## read in gene flanking gc content
with open(input_flanking_gc, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        # skip first row
        if "usercol" in row[0]:
            continue
        gene_name = row[3]
        seq_length = int(row[12])
        n_count = int(row[10])
        n_percent = n_count/seq_length
        if seq_length >= 100 and n_percent < 0.1: # make sure at least 10% of seq is actually there and make sure less than 10% of seq is Ns
            flanking_gc_count = int(row[7]) + int(row[8])
            flanking_at_count = int(row[6]) + int(row[9])
            flanking_gc_perc = flanking_gc_count/(flanking_gc_count + flanking_at_count)
            if len(gene_dict[gene_name]) == 1: # check if only one flanking region counted so far
                gene_dict[gene_name].append(flanking_gc_perc)
            elif len(gene_dict[gene_name]) == 2: # if the gene already has a flanking content number then average the two
                previous_value = gene_dict[gene_name][1]
                # calc average
                average_flanking_gc_perc = (previous_value + flanking_gc_perc)/2
                gene_dict[gene_name][1] = average_flanking_gc_perc # replace value

# filter out genes without flanking gc content
okay_genes = []
for gene in gene_dict:
    if len(gene_dict[gene]) == 2:
        okay_genes.append(gene)

with open('gc_content_table/' + genome_name + '.gc_table.txt', 'w', newline = '') as output_csv:
    w = csv.writer(output_csv, delimiter = '\t')
    for gene in gene_dict:
        if gene in okay_genes:
            w.writerow([gene, gene_dict[gene][0], gene_dict[gene][1]])