import csv
import numpy as np
import sys

pacbio_file = sys.argv[1]
illumina_file = sys.argv[2]
illumina_file_split_reads = sys.argv[3]
output = sys.argv[4]
sample = sys.argv[5]

pacbio_eccs = []
with open(pacbio_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        pacbio_eccs.append([row[0], int(row[1]), int(row[2])])

illumina = []
with open(illumina_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        illumina.append([row[0], row[1], row[2]])
illumina_indexed = [[] for i in range(56)]
for ecc in illumina:
    scaffold_num = int(ecc[0][10:12])-1
    illumina_indexed[scaffold_num].append(ecc)
illumina_arrays = []
for i in range(len(illumina_indexed)):
    illumina_arrays.append(np.array(illumina_indexed[i], dtype=object))

pacbio_with_overlap = []
overlap_count = []
pacbio_no_overlap = []
tolerance = 10
for ecc in pacbio_eccs:
    start_ecc = ecc[1]
    end_ecc = ecc[2]
    illumina_for_scaffold = illumina_arrays[int(ecc[0][10:12])-1]
    ecc_matches = illumina_for_scaffold[np.logical_and(np.isclose((illumina_for_scaffold[:,1]).astype(int), start_ecc, atol=tolerance, rtol=0),
                                    np.isclose((illumina_for_scaffold[:,2]).astype(int), end_ecc, atol=tolerance, rtol=0))]
    if np.shape(ecc_matches)[0] > 0:
        pacbio_with_overlap.append(ecc)
        overlap_count.append(np.shape(ecc_matches)[0])
    else:
        pacbio_no_overlap.append(ecc)

illumina_splitreads = []
with open(illumina_file_split_reads, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        illumina_splitreads.append([row[0], row[1], row[2]])
illumina_splitreads_indexed = [[] for i in range(56)]
for ecc in illumina_splitreads:
    scaffold_num = int(ecc[0][10:12])-1
    illumina_splitreads_indexed[scaffold_num].append(ecc)
illumina_splitreads_arrays = []
for i in range(len(illumina_splitreads_indexed)):
    illumina_splitreads_arrays.append(np.array(illumina_splitreads_indexed[i], dtype=object))

pacbio_with_overlap_splitreads = []
overlap_count = []
pacbio_no_overlap_splitreads = []
tolerance = 10
for ecc in pacbio_eccs:
    start_ecc = ecc[1]
    end_ecc = ecc[2]
    illumina_for_scaffold = illumina_splitreads_arrays[int(ecc[0][10:12])-1]
    ecc_matches = illumina_for_scaffold[np.logical_and(np.isclose((illumina_for_scaffold[:,1]).astype(int), start_ecc, atol=tolerance, rtol=0),
                                    np.isclose((illumina_for_scaffold[:,2]).astype(int), end_ecc, atol=tolerance, rtol=0))]
    if np.shape(ecc_matches)[0] > 0:
        pacbio_with_overlap_splitreads.append(ecc)
        overlap_count.append(np.shape(ecc_matches)[0])
    else:
        pacbio_no_overlap_splitreads.append(ecc)



pacbio_eccs_count = len(pacbio_eccs)
pacbio_eccs_overlap = len(pacbio_with_overlap)
pacbio_eccs_splitreads_overlap = len(pacbio_with_overlap_splitreads)

with open(output, 'w', newline = '') as output_csv:
    w = csv.writer(output_csv, delimiter = '\t')
    w.writerow([sample, pacbio_eccs_count, pacbio_eccs_overlap, pacbio_eccs_splitreads_overlap])