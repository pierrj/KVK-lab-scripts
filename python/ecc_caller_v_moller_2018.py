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

moller_2018 = sys.argv[1]
sample_moller = sys.argv[2]
ecc_caller = sys.argv[3]
ecc_caller_uniq = sys.argv[4]
tissue = sys.argv[5]
sample = sys.argv[6]
output = sys.argv[7]

## USAGE ##
# this script compares eccDNA forming regions called by Moller et al 2018 and those called by ecc_caller with the same sequencing data
# options:
# "moller_2018" - eccdnas from moller et al
# "sample_moller" - sample name in moller et al
# "ecc_caller" - ecc_caller out for sample
# "ecc_caller_uniq" - ecc_caller out for sample but only with uniquely mapped reads
# "tissue" - moller et al tissue
# "sample" - sample name
# "output" - output file

# read in moller called eccdna forming regions
moller_2018_eccs = []
with open(moller_2018, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        if row[0] == sample_moller and row[1] != 'chrMT' and row[9] != "lowqual" : ## remove low quality eccdnas
            if row[1][3:] == 'X': ## remove x and y chroms and mitochondria
                chrom = 23
            elif row[1][3:] == 'Y':
                chrom = 24
            else:
                chrom = int(row[1][3:])
            moller_2018_eccs.append([chrom, int(row[3]), int(row[4]), row[9]])

ecc_caller_eccs = []
with open(ecc_caller, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        chrom = int(row[0][7:9].strip("0"))
        if row[4] != 'lowq': ## remove lowq eccdnas
            ecc_caller_eccs.append([chrom, int(row[1]), int(row[2]), row[4]])

# index for speeding up comparison with numpy
moller_indexed = [[] for i in range(24)]
for ecc in moller_2018_eccs:
    scaffold_num = ecc[0]-1
    moller_indexed[scaffold_num].append(ecc)
moller_arrays = []
for i in range(len(moller_indexed)):
    moller_arrays.append(np.array(moller_indexed[i], dtype=object))

# count overlap
ecc_caller_eccs_with_overlap = []
ecc_caller_eccs_no_overlap = []
tolerance = 10 # overlap counts if they are within 10 bp
for ecc in ecc_caller_eccs:
    start_region = ecc[1]
    end_region= ecc[2]
    eccs_for_scaffold = moller_arrays[ecc[0]-1]
    ecc_matches = eccs_for_scaffold[np.logical_and(np.isclose((eccs_for_scaffold[:,1]).astype(int), start_region, atol=tolerance, rtol=0),
                                    np.isclose((eccs_for_scaffold[:,2]).astype(int), end_region, atol=tolerance, rtol=0))] # look to see if the two eccdnas mostly overlap
    if np.shape(ecc_matches)[0] > 0:
        ecc_caller_eccs_with_overlap.append(ecc)
    else:
        ecc_caller_eccs_no_overlap.append(ecc)


ecc_caller_eccs_highcoverage = []
with open(ecc_caller, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        chrom = int(row[0][7:9].strip("0"))
        if row[4] != 'lowq' and int(row[3]) > 10: # only high coverage region, more than 10 split reads
            ecc_caller_eccs_highcoverage.append([chrom, int(row[1]), int(row[2]), row[4]])

# count overlap
ecc_caller_eccs_highcoverage_with_overlap = []
ecc_caller_eccs_highcoverage_no_overlap = []
tolerance = 10 # overlap counts if they are within 10 bp
for ecc in ecc_caller_eccs_highcoverage:
    start_region = ecc[1]
    end_region= ecc[2]
    eccs_for_scaffold = moller_arrays[ecc[0]-1]
    ecc_matches = eccs_for_scaffold[np.logical_and(np.isclose((eccs_for_scaffold[:,1]).astype(int), start_region, atol=tolerance, rtol=0),
                                    np.isclose((eccs_for_scaffold[:,2]).astype(int), end_region, atol=tolerance, rtol=0))] # look to see if the two eccdnas mostly overlap
    if np.shape(ecc_matches)[0] > 0:
        ecc_caller_eccs_highcoverage_with_overlap.append(ecc)
    else:
        ecc_caller_eccs_highcoverage_no_overlap.append(ecc)         

# only uniq eccdnas now
ecc_caller_eccs_uniq = []
with open(ecc_caller_uniq, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        chrom = int(row[0][7:9].strip("0"))
        if row[4] != 'lowq':
            ecc_caller_eccs_uniq.append([chrom, int(row[1]), int(row[2]), row[4]])

# count overlap
ecc_caller_eccs_uniq_with_overlap = []
ecc_caller_eccs_uniq_no_overlap = []
tolerance = 10 # overlap counts if they are within 10 bp
for ecc in ecc_caller_eccs_uniq:
    start_region = ecc[1]
    end_region= ecc[2]
    eccs_for_scaffold = moller_arrays[ecc[0]-1]
    ecc_matches = eccs_for_scaffold[np.logical_and(np.isclose((eccs_for_scaffold[:,1]).astype(int), start_region, atol=tolerance, rtol=0),
                                    np.isclose((eccs_for_scaffold[:,2]).astype(int), end_region, atol=tolerance, rtol=0))] # look to see if the two eccdnas mostly overlap
    if np.shape(ecc_matches)[0] > 0:
        ecc_caller_eccs_uniq_with_overlap.append(ecc)
    else:
        ecc_caller_eccs_uniq_no_overlap.append(ecc)

## count out overlaps and output
moller_2018_eccs_count = len(moller_2018_eccs)
ecc_caller_eccs_count = len(ecc_caller_eccs)
ecc_caller_eccs_w_overlap_count = len(ecc_caller_eccs_with_overlap)
ecc_caller_eccs_highcoverage_count = len(ecc_caller_eccs_highcoverage)
ecc_caller_eccs_highcoverage_w_overlap_count = len(ecc_caller_eccs_highcoverage_with_overlap)
ecc_caller_uniq_eccs_count = len(ecc_caller_eccs_uniq)
ecc_caller_uniq_eccs_w_overlap_count = len(ecc_caller_eccs_uniq_with_overlap)

with open(output, 'w', newline = '') as output_csv:
    w = csv.writer(output_csv, delimiter = '\t')
    w.writerow([sample, tissue, moller_2018_eccs_count, ecc_caller_eccs_count, ecc_caller_eccs_w_overlap_count, 
            ecc_caller_eccs_highcoverage_count, ecc_caller_eccs_highcoverage_w_overlap_count,
            ecc_caller_uniq_eccs_count, ecc_caller_uniq_eccs_w_overlap_count])