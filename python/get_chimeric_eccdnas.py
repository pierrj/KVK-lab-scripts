import csv
import numpy as np
import sys
import re

pacbio_alignments = sys.argv[1]
illumina_splitreads = sys.argv[2]
output = sys.argv[3]
output_name = sys.argv[4]
tolerance = sys.argv[5]
column_cutoff = sys.argv[6]

def process_cigar(cigar, sense):
    start_pattern = "^([0-9]+)[HS].*[MDIHS]$"
    end_pattern = ".*[MDIHS]([0-9]+)[HS]$"
    ## need to add some plus strand and reverse strand here
    ## see m54213_190513_082505/69206095/ccs for example
    if sense == '-':
        if re.match(start_pattern, cigar):
            return int(re.match(start_pattern, cigar).group(1))
        else:
            return 0
    elif sense == '+':
        if re.match(end_pattern, cigar):
            return int(re.match(end_pattern, cigar).group(1))
        else:
            return 0

def count_matches(cigar):
    matches = re.findall(r'(\d+)([A-Z]{1})', cigar)
    matches_sums = {'M': 0, 'other': 0}
    for i in range(len(matches)):
        if matches[i][1] == 'M' or matches[i][1] == 'I':
            matches_sums['M'] += int(matches[i][0])            
        else:
            matches_sums['other'] += int(matches[i][0])
    return matches_sums['M']

reads = {}
with open(pacbio_alignments, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        if row[3] not in reads:
            reads[row[3]] = []
        reads[row[3]].append([row[0], int(row[1]), int(row[2]), int(row[4]), row[5], row[6], int(row[7]), process_cigar(row[6], row[5]), count_matches(row[6])])
reads_arrays = {}
for key in reads.keys():
    reads_arrays[key] = np.array(reads[key], dtype=object)

## filter out repetitives
unique_reads = {}
for key in reads_arrays:
    array = reads_arrays[key]
    ## sam flag 256
    if not np.array(np.bitwise_and(array[:,6], 0x100), dtype=bool).any():
        unique_reads[key] = array

chimeras = {}
for key in unique_reads:
    if np.shape(unique_reads[key])[0] > 1 and not np.all(unique_reads[key][:,0] == unique_reads[key][0][0]):
        chimeras[key] = unique_reads[key]

illumina_junctions = []
with open(illumina_splitreads, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        illumina_junctions.append([row[0], int(row[1]), row[2], int(row[3])])
illumina_junctions_arrays = np.array(illumina_junctions, dtype=object)

qced_chimeras = {}
for key in chimeras:
    chimera = chimeras[key]
    chimera_sorted = chimera[chimera[:, 7].argsort()]
    proper_order_test = 0
    for i in range(len(chimera_sorted)):
        if i != 0:
            soft_clipped = chimera_sorted[i][7]
            predicted_soft_clipped = chimera_sorted[i-1][7] + chimera_sorted[i-1][8]
            if not np.isclose(soft_clipped, predicted_soft_clipped, atol=100, rtol=0):
                proper_order_test += 1
    if proper_order_test == 0:
        scaffold_change_test = 0
        previous_alignment_scaffold = 0
        for alignment in chimera_sorted:
            if previous_alignment_scaffold == alignment[0]:
                scaffold_change_test += 1
            previous_alignment_scaffold = alignment[0]
        if scaffold_change_test == 0:
            qced_chimeras[key] = {'junctions': []}
            qced_chimeras[key]['alignments'] = chimera_sorted
            for i in range(len(chimera_sorted)):
                if i != 0:
                    if chimera_sorted[i-1][4] == '+':
                        junction_start = chimera_sorted[i-1][1]
                    elif chimera_sorted[i-1][4] == '-':
                        junction_start = chimera_sorted[i-1][2]
                    if chimera_sorted[i][4] == '+':
                        junction_end = chimera_sorted[i][2]
                    elif chimera_sorted[i][4] == '-':
                        junction_end = chimera_sorted[i][1]
                    qced_chimeras[key]['junctions'].append([ chimera_sorted[i-1][0], junction_start ,chimera_sorted[i][0], junction_end])

confirmed_chimeras = []
for key in qced_chimeras:
    junctions = qced_chimeras[key]['junctions']
    junction_check = []
    for junction in junctions:
        confirming_illumina = illumina_junctions_arrays[np.logical_or(
            np.logical_and(illumina_junctions_arrays[:,0] == junction[0],
                           illumina_junctions_arrays[:,2] == junction[2]),
            np.logical_and(illumina_junctions_arrays[:,2] == junction[0],
                           illumina_junctions_arrays[:,0] == junction[2]))]
        confirming_illumina = confirming_illumina[np.logical_or(
        np.logical_and(np.isclose((confirming_illumina[:,1]).astype(int), junction[1], atol=10, rtol=0), 
                       np.isclose((confirming_illumina[:,3]).astype(int), junction[3], atol=10, rtol=0)),
        np.logical_and(np.isclose((confirming_illumina[:,3]).astype(int), junction[1], atol=10, rtol=0),
                       np.isclose((confirming_illumina[:,1]).astype(int), junction[3], atol=10, rtol=0)))]
        if np.shape(confirming_illumina)[0] == 0:
            junction_check.append(0)
    if 0 not in junction_check:
        ## one last filter here, all alignments should be same orientation
        chimera = qced_chimeras[key]['alignments']
        all_same_orientation = 0
        for scaffold in np.unique(chimera[:,0]):
            scaffold_subset = chimera[chimera[:,0] == scaffold]
            if not np.all(scaffold_subset[:,4] == scaffold_subset[0][4]):
                all_same_orientation += 1
        if all_same_orientation == 0:
            confirmed_chimeras.append(chimera)

chimera_regions = []
for i in range(len(confirmed_chimeras)):
    chimera = confirmed_chimeras[i]
    for scaffold in np.unique(chimera[:,0]):
        scaffold_subset = chimera[chimera[:,0] == scaffold]
        min_scaffold = np.min(scaffold_subset[:,1])
        max_scaffold = np.max(scaffold_subset[:,2])
        chimera_regions.append([scaffold, min_scaffold, max_scaffold, output_name + '_confirmedchimera_'+str(i)])

with open(output, 'w', newline="") as output_file:
    w = csv.writer(output_file, delimiter = '\t')
    for row in chimera_regions:
        w.writerow([row[0], row[1], row[2], row[3]])