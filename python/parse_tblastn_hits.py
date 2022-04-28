import numpy as np
import csv
import more_itertools as mit
from itertools import filterfalse
import sys

input_file = sys.argv[1]
e_value = float(sys.argv[2])
pident = int(sys.argv[3])
query_cov = int(sys.argv[4])
hit_count = int(sys.argv[5])
genome = sys.argv[6]
og = sys.argv[7]

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
                if protein[:-2] not in valid_hits: # same protein cant be counted twice for two alignments
                    valid_hits.append(protein[:-2])

if len(protein_hits) >= hit_count:
    print(genome + '\t' + og + '\tyes')
    bed_entries = []
    for hit in valid_hits:
        scaffold = valid_hits[hit][:,0][0]
        start = np.min(valid_hits[hit][:,5:7])
        end = np.max(valid_hits[hit][:,5:7])
        bed = np.array([scaffold, start, end])
        if bed not in bed_entries:
            bed_entries.append(bed)
    for count, bed in enumerate(bed_entries):
        with open(genome+'_'+og+'_'+count+'.bed', 'w', newline = '') as output_bed:
                  w = csv.writer(output_bed,delimiter = '\t')
                  w.writewrow(bed)
else:
    print(genome + '\t' + og + '\tno')