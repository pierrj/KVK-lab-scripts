import csv
import sys
from collections import Counter

input_file = sys.argv[1]
expected_og = sys.argv[2]
genome = sys.argv[3]

hits_dict = {}

with open(input_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        if row[0] not in hits_dict:
            hits_dict[row[0]] = []
        hits_dict[row[0]].append(row[1][-9:])

weighed_counts = {}

for protein in hits_dict:
    number_of_overlapping_proteins = int(protein.split('_')[1])
    c = Counter(hits_dict['0_122'])
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

if max(weighed_counts_percents.values()) < 1.0:
    print('some disagreement in ogs, max was')
    print(print(max(weighed_counts_percents.values())))
observed_og = max(weighed_counts_percents, key=weighed_counts_percents.get)
if expected_og == observed_og:
    print(genome + '\t' + expected_og + '\tyes')
else:
    print(genome + '\t' + expected_og + '\tno_wrong_blastp_hit')