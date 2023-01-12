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
from collections import Counter

input_file = sys.argv[1]
expected_og = sys.argv[2]
genome = sys.argv[3]

expected_og = expected_og.split('.')[0]

## read in blast p hits
hits_dict = {}

with open(input_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        if row[0] not in hits_dict:
            hits_dict[row[0]] = []
        hits_dict[row[0]].append(row[1][-9:])

if not hits_dict: # sometimes there aren't any blastp hits
    print(genome + '\t' + expected_og + '\tno_blastp_hit')
    exit()

weighed_counts = {}

for protein in hits_dict:
    number_of_overlapping_proteins = int(protein.split('_')[1]) ## blastp hits, should be 100 in total
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

## pick the og with the most blastp hits
observed_og = max(weighed_counts_percents, key=weighed_counts_percents.get)

if expected_og == observed_og:
    print(genome + '\t' + expected_og + '\tyes'+'\t'+str(max(weighed_counts_percents.values()))+'\t'+observed_og)
else:
    print(genome + '\t' + expected_og + '\tno_wrong_blastp_hit'+'\t'+str(max(weighed_counts_percents.values()))+'\t'+observed_og)