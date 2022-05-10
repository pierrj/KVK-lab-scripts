import csv
import sys
from collections import Counter

input_file = sys.argv[1]
expected_og = sys.argv[2]
genome = sys.argv[3]

expected_og = expected_og.split('.')[0]

hits_dict = {}

with open(input_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        if row[0] not in hits_dict:
            hits_dict[row[0]] = []
        hits_dict[row[0]].append(row[1][-9:])

if not hits_dict: # sometimes the blastp hit is emtpy
    print(genome + '\t' + expected_og + '\tno_blastp_hit')
    exit()

weighed_counts = {}

for protein in hits_dict:
    number_of_overlapping_proteins = int(protein.split('_')[1]) ## tblastn hits
    c = Counter(hits_dict[protein])
    if len(hits_dict[protein]) < 100:
        for og in c:
            c[og] = round(c[og]*100/len(hits_dict[protein])) ## if less than 100 hits, weigh them to 100
    for og in c:
        if og not in weighed_counts:
            weighed_counts[og] = c[og]*number_of_overlapping_proteins
        else:
            weighed_counts[og] += c[og]*number_of_overlapping_proteins

weighed_counts_sum = 0

for weighed_count in weighed_counts:
    weighed_counts_sum += weighed_counts[weighed_count]

weighed_counts_percents = {}

for weighed_count in weighed_counts:
    weighed_counts_percents[weighed_count] = weighed_counts[weighed_count]/weighed_counts_sum

observed_og = max(weighed_counts_percents, key=weighed_counts_percents.get)

if expected_og == observed_og:
    print(genome + '\t' + expected_og + '\tyes'+'\t'+str(max(weighed_counts_percents.values()))+'\t'+observed_og)
else:
    print(genome + '\t' + expected_og + '\tno_wrong_blastp_hit'+'\t'+str(max(weighed_counts_percents.values()))+'\t'+observed_og)