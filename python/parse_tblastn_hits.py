import numpy as np
import csv
import more_itertools as mit
from itertools import filterfalse
import sys

input_file = sys.argv[1]
e_value = float(sys.argv[2])
pident = int(sys.argv[3])
query_cov = int(sys.argv[4])
hit_count = int(sys.argv[5])
og = sys.argv[6]

prelim_hits = {}

with open(input_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        if row[0] not in prelim_hits:
            prelim_hits[row[0]] = []
        prelim_hits[row[0]].append([row[1],float(row[2]),int(row[3]),int(row[4]),
                             int(row[5]),int(row[6]),int(row[7]),float(row[8])])
prelim_hits_arrays = {}

for key in prelim_hits.keys():
    prelim_hits_arrays[key] = np.array(prelim_hits[key], dtype=object)

parsed_hits = {}

for protein in prelim_hits.keys():
    prelim_hit_sorted = prelim_hits_arrays[protein][prelim_hits_arrays[protein][:, 5].argsort()]
    count = 0
    previous = 0
    for index in range(len(prelim_hit_sorted)):
        hsp = prelim_hit_sorted[index]
        if np.any(previous != 0):
            if hsp[5] - previous[6] < 3000:
                parsed_hits[protein+'_'+str(count)].append(hsp)
                previous = hsp
            elif index == len(prelim_hit_sorted)-1:
                count += 1
                parsed_hits[protein+'_'+str(count)] = []
                parsed_hits[protein+'_'+str(count)].append(hsp)
            else:
                count += 1
                parsed_hits[protein+'_'+str(count)] = []
                parsed_hits[protein+'_'+str(count)].append(hsp)
                previous = hsp
        else:
            previous = hsp
            parsed_hits[protein+'_'+str(count)] = []
            parsed_hits[protein+'_'+str(count)].append(hsp)
            
parsed_hits_arrays = {}

for key in parsed_hits.keys():
    parsed_hits_arrays[key] = np.array(parsed_hits[key], dtype=object)

valid_hits = []

for protein in parsed_hits_arrays:
    hit = parsed_hits_arrays[protein]
    if (np.max(hit[:,0]) == np.min(hit[:,0]) and # all same scaffold
       np.all(hit[:,1] < e_value) and # e-value cutoff
       np.all(hit[:,7] > pident) # pident cutoff
       ):
        protein_size_range = range(1,hit[0,2])
        for i in hit:
            protein_size_range = list(filterfalse(lambda x: i[3] <= x <= i[4], protein_size_range)) # get query cov
        if (1-len(protein_size_range))*100 > query_cov: # check if query cov is enough
            if protein[:-2] not in valid_hits: # same protein cant be counted twice for two alignments
                valid_hits.append(protein)

if len(valid_hits) >= hit_count:
    print(og + '\tyes')
else:
    print(og+'\tno')

same_scaffold_filtered_hits = []
for protein in parsed_hits_arrays:
    hit = parsed_hits_arrays[protein]
    if np.max(hit[:,0]) == np.min(hit[:,0]):
        if protein[:-2] not in same_scaffold_filtered_hits:
            same_scaffold_filtered_hits.append(protein[:-2])

evalue_filtered_hits = []
for protein in parsed_hits_arrays:
    hit = parsed_hits_arrays[protein]
    if np.all(hit[:,1] < e_value):
        if protein[:-2] not in evalue_filtered_hits:
            evalue_filtered_hits.append(protein[:-2])

pident_filtered_hits = []
for protein in parsed_hits_arrays:
    hit = parsed_hits_arrays[protein]
    if np.all(hit[:,7] > pident):
        if protein[:-2] not in pident_filtered_hits:
            pident_filtered_hits.append(protein[:-2])

query_cov_filtered_hits = []
for protein in parsed_hits_arrays:
    protein_size = hit[0,2]
    hit = parsed_hits_arrays[protein]
    protein_size_range = range(1,protein_size)
    for i in hit:
        protein_size_range = list(filterfalse(lambda x: i[3] <= x <= i[4], protein_size_range)) # get query cov
    if (1-(len(protein_size_range)/protein_size))*100 > query_cov:
        if protein[:-2] not in query_cov_filtered_hits:
            query_cov_filtered_hits.append(protein[:-2])

print('total')
print(len(list(parsed_hits_arrays.keys())))
print('same scaffold filter')
print(len(same_scaffold_filtered_hits))
print('evalue filter')
print(len(evalue_filtered_hits))
print('pident filter')
print(len(pident_filtered_hits))
print('query cov filter')
print(len(query_cov_filtered_hits))