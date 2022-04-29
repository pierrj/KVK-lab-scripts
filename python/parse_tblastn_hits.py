import numpy as np
import csv
import more_itertools as mit
from itertools import filterfalse
import sys
import os

input_file = sys.argv[1]
e_value = float(sys.argv[2])
pident = int(sys.argv[3])
query_cov = int(sys.argv[4])
hit_count = int(sys.argv[5])
genome = sys.argv[6]
og = sys.argv[7]
gene_gff = sys.argv[8]
gene_overlap = float(sys.argv[9])

prelim_hits = {}

with open(input_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        if row[0] not in prelim_hits:
            prelim_hits[row[0]] = []
        prelim_hits[row[0]].append([row[1],float(row[2]),int(row[3]),int(row[4]),
                             int(row[5]),int(row[6]),int(row[7]),
                             (int(row[8])/(int(row[8])+int(row[9])))*100]) # calculate pident manually here
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

protein_hits = []
valid_hits = []

for protein in parsed_hits_arrays:
    hit = parsed_hits_arrays[protein]
    hit = hit[hit[:,1] < e_value] # remove hsps below evalue
    hit = hit[hit[:,7] > pident] # remove hsps below pident
    if hit.size != 0: # check that it isn't empty
        if np.max(hit[:,0]) == np.min(hit[:,0]): # make sure all are from same scaffold
            protein_size = hit[0,2]
            protein_size_range = range(1,protein_size)
            for i in hit:
                protein_size_range = list(filterfalse(lambda x: i[3] <= x <= i[4], protein_size_range)) # get query cov
            if (1-(len(protein_size_range)/protein_size))*100 > query_cov: # check if query cov for remaining hsps is enough
                valid_hits.append(hit)
                if protein[:-2] not in protein_hits: # same protein cant be counted twice for two alignments
                    protein_hits.append(protein[:-2])


def numpy_get_overlap_percent(lst, twod_array, percent): ## function to quickly calculate percent overlap between two alignments from the input coords file
    twod_array_min = twod_array.min(axis=1)
    twod_array_max = twod_array.max(axis=1)
    twod_array_length = twod_array_max-twod_array_min
    twod_array_ordered = np.vstack((twod_array_min, twod_array_max)).T # transpose
    lst_start = min(lst)
    lst_end = max(lst)
    number_of_matches = np.shape(twod_array)[0] # count matches
    lst_start_array = np.repeat(lst_start, number_of_matches)
    lst_end_array = np.repeat(lst_end, number_of_matches)
    lst_length = lst_end-lst_start
    zero_array = np.repeat(0, number_of_matches) # count zeroes
    ## watch out for integers here
    overlap_array = np.vstack((lst_end_array, twod_array_ordered[:,1])).T.min(axis=1) - np.vstack((lst_start_array, twod_array_ordered[:,0])).T.max(axis=1)
    overlap_or_zero = np.vstack((zero_array, overlap_array)).T.max(axis=1)
    overlap_percent = overlap_or_zero/lst_length
    print(overlap_percent)
    boolean_array = overlap_percent > percent # compare to percentage cutoff
    return boolean_array

if len(protein_hits) >= hit_count:
    gene_gff_list = []
    with open(gene_gff, newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        for row in file_reader:
            gene_gff_list.append([row[0], int(row[3]), int(row[4])])
    gene_gff_arrays = np.array(gene_gff_list)
    bed_entries = []
    for hit in valid_hits:
        scaffold = hit[:,0][0]
        start = np.min(hit[:,5:7])
        end = np.max(hit[:,5:7])
        bed = [scaffold, start, end]
        if bed not in bed_entries:
            bed_entries.append(bed)
    beds_with_intersect = []
    for bed in bed_entries:
        print(bed)
        subset = gene_gff_arrays[gene_gff_arrays[:,0] == bed[0]]
        subset = subset[
            np.logical_or(
                np.logical_and(
                    subset[:,1].astype(int) <= bed[1],
                    subset[:,2].astype(int) >= bed[1]
                ),
                np.logical_and(
                    subset[:,1].astype(int) <= bed[2],
                    subset[:,2].astype(int) >= bed[2]
                )
            )
        ]
        subset = numpy_get_overlap_percent([bed[1], bed[2]],subset[:,1:3].astype(int), gene_overlap)
        if subset.size != 0:
            beds_with_intersect.append(bed)
    if len(beds_with_intersect) > 0:
            print(genome + '\t' + og + '\tyes_otherog')
    else:
        print(genome + '\t' + og + '\tyes_unannotated')
else:
    print(genome + '\t' + og + '\tno')