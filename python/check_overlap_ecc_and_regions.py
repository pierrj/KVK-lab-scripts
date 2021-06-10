import csv
import numpy as np
import sys

region_file = sys.argv[1]
ecc_file = sys.argv[2]
tolerance = int(sys.argv[3])

regions = []
with open(region_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        regions.append([row[0], int(row[1]), int(row[2])])

eccs = []
with open(ecc_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        eccs.append([row[0], row[1], row[2]])


eccs_indexed = [[] for i in range(56)]
for ecc in eccs:
    scaffold_num = int(ecc[0][10:12])-1
    eccs_indexed[scaffold_num].append(ecc)
eccs_arrays = []
for i in range(len(eccs_indexed)):
    eccs_arrays.append(np.array(eccs_indexed[i], dtype=object))

regions_with_overlap = []
for region in regions:
    start_region = region[1]
    end_region= region[2]
    eccs_for_scaffold = eccs_arrays[int(region[0][10:12])-1]
    ecc_matches = eccs_for_scaffold[np.logical_and(np.isclose((eccs_for_scaffold[:,1]).astype(int), start_region, atol=tolerance, rtol=0),
                                    np.isclose((eccs_for_scaffold[:,2]).astype(int), end_region, atol=tolerance, rtol=0))]
    if np.shape(ecc_matches)[0] > 0:
        regions_with_overlap.append(ecc_matches)

print(len(regions_with_overlap))