# coding: utf-8

# In[1]:


import csv
import ipyparallel as ipp
from itertools import groupby
from itertools import compress
import subprocess
import statistics
import collections
import sys
import os
from operator import itemgetter

with open('/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/results/ecccaller_output.G3_1A.details.tsv', newline = '') as eccloc:
    eccloc_reader = csv.reader(eccloc, delimiter = '\t')
    G3_1A_eccs = [[int(row[0][10:12]) - 1, int(row[1]), int(row[2]), int(row[3]), 3] for row in eccloc_reader]

with open('/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/results/ecccaller_output.G3_1B.details.tsv', newline = '') as eccloc:
    eccloc_reader = csv.reader(eccloc, delimiter = '\t')
    G3_1B_eccs = [[int(row[0][10:12]) - 1, int(row[1]), int(row[2]), int(row[3]), 4] for row in eccloc_reader]


with open('/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/results/ecccaller_output.G3_1C.details.tsv', newline = '') as eccloc:
    eccloc_reader = csv.reader(eccloc, delimiter = '\t')
    G3_1C_eccs = [[int(row[0][10:12]) - 1, int(row[1]), int(row[2]), int(row[3]), 5] for row in eccloc_reader]

G3_1_eccs = G3_1A_eccs + G3_1B_eccs + G3_1C_eccs

len(G3_1_eccs)

G3_1_eccs_sorted = sorted(G3_1_eccs, key=itemgetter(0,1))

print(G3_1_eccs_sorted[0:50])

chrom10 = []
for i in range(len(G3_1_eccs_sorted)):
    if G3_1_eccs_sorted[i][0] == 9:
        chrom10.append(G3_1_eccs_sorted[i])

loc_list = [[] for i in range(56)]
target = G3_1_eccs_sorted
for i in range(len(target)):
    loc = [target[i][0], target[i][1], target[i][2]]
    loc_list[target[i][0]].append(loc)

loc_list_uniq = [sorted([list(x) for x in set(tuple(x) for x in loc_list[k])]) for k in range(len(loc_list))]

for i in range(len(loc_list_uniq)):
    for l in range(len(loc_list_uniq[i])):
        loc_list_uniq[i][l].append(0) ## G3_1A index is 3
        loc_list_uniq[i][l].append(0) ## G3_1B index is 4
        loc_list_uniq[i][l].append(0) ## G3_1C index is 5

for i in range(len(target)):
    ecc = target[i]
    for l in range(len(loc_list_uniq[ecc[0]])):
        if ecc[1] == loc_list_uniq[ecc[0]][l][1] and ecc[2] == loc_list_uniq[ecc[0]][l][2]:
            loc_list_uniq[ecc[0]][l][ecc[4]] += ecc[3]

for l in range(len(loc_list_uniq)):
    for i in range(len(loc_list_uniq[l])):
        if loc_list_uniq[l][i][3] > 0 and loc_list_uniq[l][i][4] > 0 and loc_list_uniq[l][i][5] > 0:
            print(loc_list_uniq[l][i])

def get_variants_grouped(eccs_perchrom):
    ecc_withvariants = {}
    for i in range(len(eccs_perchrom)):
        ecc = eccs_perchrom[i]
        variants = []
# rolling grouping of eccs by proximity
        for k in range(len(eccs_perchrom)):
            start_coordinate1 = eccs_perchrom[k][1] - 20
            start_coordinate2 = eccs_perchrom[k][1] + 20
            end_coordinate1 = eccs_perchrom[k][2] - 20
            end_coordinate2 = eccs_perchrom[k][2] + 20
            if start_coordinate1 <= ecc[1] <= start_coordinate2 and end_coordinate1 <= ecc[2] <= end_coordinate2:
                variants.append(eccs_perchrom[k])
        # check if there isnt just itself in the list
        if tuple([ecc]) != tuple(variants):
            ecc_withvariants[tuple(ecc)] = variants
    variants_grouped = collections.defaultdict(list)
    # group eccs with exact same variants
    # #because of rolling grouping this results in some less than ideal results when two groups of eccDNAs are near each other, but it only makes the grouping less accurate and isnt particularly dangerous
    for key,val in ecc_withvariants.items():
        variants_grouped[tuple(tuple(x) for x in val)].append(key)
    return variants_grouped

def getcoord(coord_list):
    median = statistics.median(coord_list)
    distancetoend = abs(median - max(coord_list))
    distancetostart = abs(median - min(coord_list))
    greaterlist = [i for i in coord_list if i >= median]
    smallerlist = [i for i in coord_list if i <= median]
    # if median is closer to end, get value above median, if median is closer to start, get value below median
    # if equal get both values, one is set to alt
    if distancetoend > distancetostart:
        coord = min(smallerlist, key=lambda x:abs(x-median))
        coord_alt = 'N/A'
    elif distancetoend == distancetostart:
        coord = min(greaterlist, key=lambda x:abs(x-median))
        coord_alt = min(smallerlist, key=lambda x:abs(x-median))
    elif distancetoend < distancetostart:
        coord = min(greaterlist, key=lambda x:abs(x-median))
        coord_alt = 'N/A'
    else:
        raise ValueError("getcoord error")
    return coord, coord_alt

def test_coords(target_coord, other_coord, target_start_or_end, val):
        # check if start or end coordinate is the one that is set by the function
        if target_start_or_end == 'start':
            known = 1
            unknown = 2
        if target_start_or_end == 'end':
            known = 2
            unknown = 1
        # get all potential unset coordinates that match set coordinate
        unknown_options = []
        for i in val:
            if i[known] == target_coord:
                unknown_options.append(i[unknown])
        # test which of the potential unset coordinates is closest to the ideal coordinate
        unknown_coord = min(unknown_options, key=lambda x:abs(x-other_coord))
        return unknown_coord


def reconcile_coords(startcoord, endcoord, val):
        # test with start coordinate set/forced
        startcoord_forced = test_coords(startcoord, endcoord, 'start', val)
        # test with end coordinate set/forced
        endcoord_forced = test_coords(endcoord, startcoord, 'end', val)
        # test how close each non-ideal corresponding coordinates are to the ideal
        fit_startcoord_forced = abs(startcoord_forced - endcoord)
        fit_endcoord_forced = abs(endcoord_forced - startcoord)
        if fit_startcoord_forced <= fit_endcoord_forced:
            return startcoord, startcoord_forced
        else:
            return endcoord_forced, endcoord

def get_representative_variants(variants_grouped):
    representative_variants = {}
    for key, val in variants_grouped.items():
        startlist = []
        endlist = []
        for i in val:
            startlist.append(i[1])
            endlist.append(i[2])
        start_coord, start_coord_alt = getcoord(startlist)
        end_coord, end_coord_alt = getcoord(endlist)
        middle_ecc = 0
    # check to see if a single ecc has ideal coordinates for both end and start, check alternates as well
        for i in val:
            if i[1] == start_coord and i[2] == end_coord:
                middle_ecc = i
            if i[1] == start_coord and i[2] == end_coord_alt:
                middle_ecc = i
            if i[1] == start_coord_alt and i[2] == end_coord:
                middle_ecc = i
            if i[1] == start_coord and i[2] == start_coord_alt:
                middle_ecc = i
    # otherwise reconcile coordinates to get closest representative ecc
        if middle_ecc == 0:
            start_coord_estimate, end_coord_estimate = reconcile_coords(start_coord, end_coord, val)
            for i in val:
                if i[1] == start_coord_estimate and i[2] == end_coord_estimate:
                    middle_ecc = i
    # make sure that a variant always gets called
        if middle_ecc == 0:
            raise ValueError("No ecc merged ecc called")
        else:
            representative_variants[middle_ecc] = val
    return representative_variants

def uniq_sort(start_list):
    uniq_list = [list(x) for x in set(tuple(x) for x in start_list)]
    end_list = sorted(uniq_list,key=lambda x: x[1])
    return end_list

def merge_variants(eccs_perchrom, representative_variants):
    for i in range(len(eccs_perchrom)):
        ecc = eccs_perchrom[i]
        for key, val in representative_variants.items():
            if tuple(ecc) in val and [tuple(ecc)] != val: #double check to make sure val isnt just ecc
                key_list = list(key)
                # replace the number of split reads by the sum of split reads of the eccDNAs
                sr_counts = [0 for count in range(3, len(val[0]))]
                for l in range(len(val)):
                    for g, k in enumerate(range(3, len(val[l]))):
                        sr_counts[g] += val[l][k]
                for x, z in enumerate(range(3, len(key_list))):
                    key_list[z] = sr_counts[x]
                eccs_perchrom[i] = key_list
    #get rid of all duplicates caused by variants being replaced by their representative eccs and then sort them
    variants_merged = uniq_sort(eccs_perchrom)
    return variants_merged

# wrapper function for merging eccDNAs into representative variants
def ecc_merge(eccs_perchrom):
    variants_grouped = get_variants_grouped(eccs_perchrom)
    representative_variants = get_representative_variants(variants_grouped)
    variants_merged = merge_variants(eccs_perchrom, representative_variants)
    return variants_merged, representative_variants

test_list = loc_list_uniq[9]
variants = get_variants_grouped(test_list)

representative_variants = get_representative_variants(variants)

variants_merged = merge_variants(test_list, representative_variants)

for i in range(len(variants_merged)):
    if variants_merged[i][3] > 0 and variants_merged[i][4] > 0 and variants_merged[i][5] > 0:
        print(variants_merged[i])


merged_list = []
splitreads_list = []
for i in range(len(loc_list_uniq)):
    variants_merged, representative_variants = ecc_merge(loc_list_uniq[i])
    merged_list.append(variants_merged)
    splitreads_list.append(representative_variants)

for l in range(len(merged_list)):
    for i in range(len(merged_list[l])):
        if merged_list[l][i][3] > 0 and merged_list[l][i][4] > 0 and merged_list[l][i][5] > 0:
            print(merged_list[l][i])

sample = 'IF_3'
with open('/global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/' + sample +'A/parallel.confirmed') as confirmed:
    confirmed_reader = csv.reader(confirmed, delimiter = '\t')
    confirmed_list_A = [[int(row[0]), int(row[1]), int(row[2])] for row in confirmed_reader]

with open('/global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/' + sample +'B/parallel.confirmed') as confirmed:
    confirmed_reader = csv.reader(confirmed, delimiter = '\t')
    confirmed_list_B = [[int(row[0]), int(row[1]), int(row[2])] for row in confirmed_reader]

with open('/global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/' + sample +'C/parallel.confirmed') as confirmed:
    confirmed_reader = csv.reader(confirmed, delimiter = '\t')
    confirmed_list_C = [[int(row[0]), int(row[1]), int(row[2])] for row in confirmed_reader]

confirmed_list = confirmed_list_A + confirmed_list_B + confirmed_list_C

confirmed_list_sorted = sorted(confirmed_list, key=itemgetter(0,1))

with open('/global/scratch/users/pierrj/eccDNA/pipeline_tests/merging/parallel.confirmed.' + sample + '.bed', 'w', newline = '') as bed:
    w = csv.writer(bed, delimiter = '\t')
    for i in range(len(confirmed_list_sorted)):
        scaffold_string = 'MQOP010000' + str(confirmed_list_sorted[i][0]+1).zfill(2) + '.1'
        row = [scaffold_string, confirmed_list_sorted[i][1], confirmed_list_sorted[i][2], i]
        w.writerow(row)



print(confirmed_list_sr_sorted[0:50])