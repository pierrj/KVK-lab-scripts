import csv
import numpy as np
import sys
import re

pacbio_alignments = sys.argv[1]
output = sys.argv[2]
tolerance = int(sys.argv[3])
column_cutoff = int(sys.argv[4])

# use regex to grab the matches and nonmatches to the genome and count them
def process_cigar(cigar, sense):
    matches = re.findall(r'(\d+)([A-Z]{1})', cigar)
    matches_sums = {'M': 0, 'other': 0}
    for i in range(len(matches)):
        if matches[i][1] == 'M':
            matches_sums['M'] += int(matches[i][0])
        else:
            matches_sums['other'] += int(matches[i][0])
    start_pattern = "^([0-9]+)[HS].*[MDIHS]$"
    end_pattern = ".*[MDIHS]([0-9]+)[HS]$"
    if sense == '+':
        if re.match(start_pattern, cigar) and re.match(end_pattern, cigar):
            if int(re.match(start_pattern, cigar).group(1)) < int(re.match(end_pattern, cigar).group(1)):
                loc = 'start'
            elif int(re.match(start_pattern, cigar).group(1)) > int(re.match(end_pattern, cigar).group(1)):
                loc = 'end'
            else:
                return False
        else:
            if re.match(start_pattern, cigar):
                loc = 'end'
            elif re.match(end_pattern, cigar):
                loc = 'start'
            else:
                return False
    elif sense == '-':
        if re.match(start_pattern, cigar) and re.match(end_pattern, cigar):
            if int(re.match(start_pattern, cigar).group(1)) < int(re.match(end_pattern, cigar).group(1)):
                loc = 'end'
            elif int(re.match(start_pattern, cigar).group(1)) > int(re.match(end_pattern, cigar).group(1)):
                loc = 'start'
            else:
                return False
        else:
            if re.match(start_pattern, cigar):
                loc = 'start'
            elif re.match(end_pattern, cigar):
                loc = 'end'
            else:
                return False
    return loc, matches_sums['M'], matches_sums['other']

reads = {}
with open(pacbio_alignments, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        if row[3] not in reads:
            reads[row[3]] = []
        reads[row[3]].append([row[0], int(row[1]), int(row[2]), int(row[4]), row[5], row[6], int(row[7])])
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

## look for single read that overlaps same spot multiple times
overlapping_unique_reads = {}
split_reads = {}
for key in unique_reads:
    alignments = unique_reads[key]
    ## check all same scaffold
    ## check all same orientation
    ## check not just one
    ## this is messy but it looks for reads that wrap around the same region
    ## either it maps within tolerance to the same region two or more times
    ## or it maps within tolerance to the same start n times and end n-1 times
    ## or it maps within tolerance to the same start n-1 times and end n-1 times
    if np.all(alignments[:,0] == alignments[0][0]) and np.all(alignments[:,4] == alignments[0][4]) and np.shape(alignments)[0] > 2:
        if ((np.count_nonzero(np.isclose((alignments[:,1]).astype(int), np.min(alignments[:,1]), atol=tolerance, rtol=0)) == np.shape(alignments)[0]-1 and
            np.count_nonzero(np.isclose((alignments[:,2]).astype(int), np.max(alignments[:,2]), atol=tolerance, rtol=0)) == np.shape(alignments)[0]-1) or
        (np.count_nonzero(np.isclose((alignments[:,1]).astype(int), np.min(alignments[:,1]), atol=tolerance, rtol=0)) == np.shape(alignments)[0] and
            np.count_nonzero(np.isclose((alignments[:,2]).astype(int), np.max(alignments[:,2]), atol=tolerance, rtol=0)) == np.shape(alignments)[0]-1) or
           (np.count_nonzero(np.isclose((alignments[:,1]).astype(int), np.min(alignments[:,1]), atol=tolerance, rtol=0)) == np.shape(alignments)[0]-1 and
            np.count_nonzero(np.isclose((alignments[:,2]).astype(int), np.max(alignments[:,2]), atol=tolerance, rtol=0)) == np.shape(alignments)[0]) or
    (np.count_nonzero(np.isclose((alignments[:,1]).astype(int), np.min(alignments[:,1]), atol=tolerance, rtol=0)) == np.shape(alignments)[0] and
            np.count_nonzero(np.isclose((alignments[:,2]).astype(int), np.max(alignments[:,2]), atol=tolerance, rtol=0)) == np.shape(alignments)[0])):
            overlapping_unique_reads[key] = alignments
    elif np.all(alignments[:,0] == alignments[0][0]) and np.all(alignments[:,4] == alignments[0][4]) and np.shape(alignments)[0] == 2:
        if alignments[0][2] > alignments[1][1]: ## if overlapping
            overlapping_unique_reads[key] = alignments
        elif abs(alignments[0][1] - alignments[1][2]) < column_cutoff: ## check if too long
            alignment_1 = alignments[0]
            alignment_2 = alignments[1]
            ## make sure i can parse the string
            alignment_1_matches = process_cigar(alignment_1[5], alignment_1[4])[1]
            alignment_2_nonmatches = process_cigar(alignment_2[5], alignment_2[4])[2]
            alignment_1_read_part = process_cigar(alignment_1[5], alignment_1[4])[0]
            alignment_2_read_part = process_cigar(alignment_2[5], alignment_2[4])[0]
            ## make sure read only matched to the two locations and nothing else
            ## and that one is called start and the other end
            if np.isclose(alignment_1_matches, alignment_2_nonmatches, atol=50, rtol=0) and alignment_1_read_part != alignment_2_read_part:
                ## check for proper split read orientation to differentiate between eccdna and intron
                if ((alignment_1_read_part == 'start' and alignment_1[4] == '+' and alignment_1[1] > alignment_2[1]) or 
                   (alignment_1_read_part == 'end' and alignment_1[4] == '-' and alignment_1[1] > alignment_2[1]) or
                   (alignment_1_read_part == 'start' and alignment_1[4] == '-' and alignment_1[1] < alignment_2[1]) or
                   (alignment_1_read_part == 'end' and alignment_1[4] == '+' and alignment_1[1] < alignment_2[1])):
                    split_reads[key] = alignments

eccdna_locs = []
for key in overlapping_unique_reads:
    alignments = overlapping_unique_reads[key]
    scaffold = alignments[0][0]
    start = np.min(alignments[:,1])
    stop = np.max(alignments[:,2])
    eccdna_locs.append([scaffold, start, stop, alignments])
for key in split_reads:
    alignments = split_reads[key]
    scaffold = alignments[0][0]
    start = np.min(alignments[:,1])
    stop = np.max(alignments[:,2])
    eccdna_locs.append([scaffold, start, stop, alignments])

with open(output, 'w', newline = '') as output_csv:
    w = csv.writer(output_csv, delimiter = '\t') ## lineterminator="\n" might fix issues with windows
    for row in eccdna_locs:
        w.writerow([row[0], row[1], row[2]])