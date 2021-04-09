import csv
import numpy as np
import sys
from statistics import median

deletions_file = sys.argv[1]
ecc_file = sys.argv[2]

deletions = []
with open(deletions_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        deletions.append([row[0], int(row[1]), int(row[2])])

eccs = []
with open(ecc_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        eccs.append([row[0], int(row[1]), int(row[2])])


eccs_indexed = [[] for i in range(56)]
for ecc in eccs:
    scaffold_num = int(ecc[0][10:12])-1
    eccs_indexed[scaffold_num].append(ecc)
eccs_arrays = []
for i in range(len(eccs_indexed)):
    eccs_arrays.append(np.array(eccs_indexed[i], dtype=object))

ecc_count_per_deletion = []
tolerance = 100
for deletion in deletions:
    start_deletion = deletion[1]
    end_deletion= deletion[2]
    eccs_for_scaffold = eccs_arrays[int(deletion[0][10:12])-1]
    ecc_matches = eccs_for_scaffold[np.logical_or(np.logical_and(eccs_for_scaffold[:,1] >= start_deletion, eccs_for_scaffold[:,1] <= end_deletion)
                                    ,np.logical_and(eccs_for_scaffold[:,2] >= start_deletion, eccs_for_scaffold[:,2] <= end_deletion))]
    ecc_count_per_deletion.append(np.shape(ecc_matches)[0])

print(median(ecc_count_per_deletion))