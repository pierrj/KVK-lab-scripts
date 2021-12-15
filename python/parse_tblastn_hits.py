import numpy as np
import csv
import sys

input_file = sys.argv[1]
e_value = float(sys.argv[2])
pident = int(sys.argv[3])
query_cov = int(sys.argv[4])
hit_count = int(sys.argv[5])
og = sys.argv[6]

hits = {}

with open(input_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        if row[0] not in hits:
            hits[row[0]] = []
        hits[row[0]].append([row[1],float(row[2]),int(row[3]),int(row[4]),
                             int(row[5]),int(row[6]),int(row[7]),float(row[8]),int(row[9]),int(row[10])])
hits_arrays = {}

for key in hits.keys():
    hits_arrays[key] = np.array(hits[key], dtype=object)

valid_hits = []

for protein in hits_arrays:
    hit = hits_arrays[protein]
    if (np.max(hit[:,0]) == np.min(hit[:,0]) and # all same scaffold
       np.all(hit[:,1] < e_value) and # e-value cutoff
       np.all(hit[:,7] > pident) and # pident cutoff
       np.all(hit[:,8] > query_cov) # query cov cutoff 
       ):
        valid_hits.append(protein)

if len(valid_hits) >= hit_count:
    print(og + '\tyes')
else:
    print(og+'\tno')