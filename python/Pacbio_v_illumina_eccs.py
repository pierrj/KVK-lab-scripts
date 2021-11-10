#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
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
import numpy as np
import sys

pacbio_file = sys.argv[1]
illumina_file = sys.argv[2]
illumina_file_split_reads = sys.argv[3]
output = sys.argv[4]
sample = sys.argv[5]

## USAGE ##
# this script takes illumina called eccdnas, illumina split reads and compares them to pacbio called eccdnas
# options:
# "pacbio_file" - eccdnas called with pacbio reads
# "illumina_file" - eccdnas called with illumina reads
# "illumina_file_split_reads" - split reads with very little qc (all false positives)
# "output" - output file
# "sample" - sample name

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