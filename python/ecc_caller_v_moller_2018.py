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



moller_2018_eccs = []
with open(moller_2018, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        if row[0] == sample_moller and row[1] != 'chrMT' and row[9] != "lowqual" :
            if row[1][3:] == 'X':
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
        if row[4] != 'lowq':
            ecc_caller_eccs.append([chrom, int(row[1]), int(row[2]), row[4]])


moller_indexed = [[] for i in range(24)]
for ecc in moller_2018_eccs:
    scaffold_num = ecc[0]-1
    moller_indexed[scaffold_num].append(ecc)
moller_arrays = []
for i in range(len(moller_indexed)):
    moller_arrays.append(np.array(moller_indexed[i], dtype=object))

ecc_caller_eccs_with_overlap = []
ecc_caller_eccs_no_overlap = []
tolerance = 10
for ecc in ecc_caller_eccs:
    start_region = ecc[1]
    end_region= ecc[2]
    eccs_for_scaffold = moller_arrays[ecc[0]-1]
    ecc_matches = eccs_for_scaffold[np.logical_and(np.isclose((eccs_for_scaffold[:,1]).astype(int), start_region, atol=tolerance, rtol=0),
                                    np.isclose((eccs_for_scaffold[:,2]).astype(int), end_region, atol=tolerance, rtol=0))]
    if np.shape(ecc_matches)[0] > 0:
        ecc_caller_eccs_with_overlap.append(ecc)
    else:
        ecc_caller_eccs_no_overlap.append(ecc)


ecc_caller_eccs_highcoverage = []
with open(ecc_caller, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        chrom = int(row[0][7:9].strip("0"))
        if row[4] != 'lowq' and int(row[3]) > 10:
            ecc_caller_eccs_highcoverage.append([chrom, int(row[1]), int(row[2]), row[4]])

ecc_caller_eccs_highcoverage_with_overlap = []
ecc_caller_eccs_highcoverage_no_overlap = []
tolerance = 10
for ecc in ecc_caller_eccs_highcoverage:
    start_region = ecc[1]
    end_region= ecc[2]
    eccs_for_scaffold = moller_arrays[ecc[0]-1]
    ecc_matches = eccs_for_scaffold[np.logical_and(np.isclose((eccs_for_scaffold[:,1]).astype(int), start_region, atol=tolerance, rtol=0),
                                    np.isclose((eccs_for_scaffold[:,2]).astype(int), end_region, atol=tolerance, rtol=0))]
    if np.shape(ecc_matches)[0] > 0:
        ecc_caller_eccs_highcoverage_with_overlap.append(ecc)
    else:
        ecc_caller_eccs_highcoverage_no_overlap.append(ecc)         

ecc_caller_eccs_uniq = []
with open(ecc_caller_uniq, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        chrom = int(row[0][7:9].strip("0"))
        if row[4] != 'lowq':
            ecc_caller_eccs_uniq.append([chrom, int(row[1]), int(row[2]), row[4]])


ecc_caller_eccs_uniq_with_overlap = []
ecc_caller_eccs_uniq_no_overlap = []
tolerance = 10
for ecc in ecc_caller_eccs_uniq:
    start_region = ecc[1]
    end_region= ecc[2]
    eccs_for_scaffold = moller_arrays[ecc[0]-1]
    ecc_matches = eccs_for_scaffold[np.logical_and(np.isclose((eccs_for_scaffold[:,1]).astype(int), start_region, atol=tolerance, rtol=0),
                                    np.isclose((eccs_for_scaffold[:,2]).astype(int), end_region, atol=tolerance, rtol=0))]
    if np.shape(ecc_matches)[0] > 0:
        ecc_caller_eccs_uniq_with_overlap.append(ecc)
    else:
        ecc_caller_eccs_uniq_no_overlap.append(ecc)

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