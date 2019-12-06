import csv
import collections
import statistics

with open('ecc_confidence') as confidence:
    confidence_reader = csv.reader(confidence, delimiter = '\t')
    confidence_list = [[] for i in range(56)]
    for row in confidence_reader:
        confidence_list[int(row[0])].append([int(row[0]), int(row[1]), int(row[2]), str(row[3]), str(row[4]), str(row[5]), 'no'])

def get_variants_grouped(eccs_perchrom):
    ecc_withvariants = {}
    for i in range(len(eccs_perchrom)):
        ecc = eccs_perchrom[i]
        variants = []
## this likely isnt ideal because of its rolling nature but this will be replaced by repeat/junction based later
        for k in range(len(eccs_perchrom)):
            start_coordinate1 = eccs_perchrom[k][1] - 50
            start_coordinate2 = eccs_perchrom[k][1] + 50
            end_coordinate1 = eccs_perchrom[k][2] - 50
            end_coordinate2 = eccs_perchrom[k][2] + 50
            if start_coordinate1 <= ecc[1] <= start_coordinate2 and end_coordinate1 <= ecc[2] <= end_coordinate2:
                variants.append(eccs_perchrom[k])
        if tuple([ecc]) != tuple(variants): ## if there isn't just itself in the list
            ecc_withvariants[tuple(ecc)] = variants
    variants_grouped = collections.defaultdict(list)
    for key,val in ecc_withvariants.items(): ## group eccs with same variants
        variants_grouped[tuple(tuple(x) for x in val)].append(key)
    return variants_grouped

def getcoord(coord_list): ## get ideal representative coordinates based off medians
    median = statistics.median(coord_list)
    distancetoend = abs(median - max(coord_list))
    distancetostart = abs(median - min(coord_list))
    greaterlist = [i for i in coord_list if i >= median]
    smallerlist = [i for i in coord_list if i <= median]
    #if median is closer to end, get value above median. if median is closer to start, get value below median
    #if equal get both values, one as alt
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

def test_coords(target_coord, other_coord, target_start_or_end, val): ## get closest matching coordinate based off on set/forced coordinate
        ## check if start or end coordinate is the set one
        if target_start_or_end == 'start':
            known = 1
            unknown = 2
        if target_start_or_end == 'end':
            known = 2
            unknown = 1
        ## get all potential options matching set coordinate
        unknown_options = []
        for i in val:
            if i[known] == target_coord:
                unknown_options.append(i[unknown])
        ## test which one is the closest to the ideal median
        unknown_coord = min(unknown_options, key=lambda x:abs(x-other_coord))
        return unknown_coord

def reconcile_coords(startcoord, endcoord, val): ## test which coordinate pair is closest to ideal
        #test with start coordinate forced
        startcoord_forced = test_coords(startcoord, endcoord, 'start', val)
        #test with end coordinate forced
        endcoord_forced = test_coords(endcoord, startcoord, 'end', val)
        #test how close each non-ideal corresponding coordinates are to the ideal
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
    ## check to see if a single ecc has ideal coordinates for both end and start
        for i in val:
            if i[1] == start_coord and i[2] == end_coord:
                middle_ecc = i
            if i[1] == start_coord and i[2] == end_coord_alt:
                middle_ecc = i
            if i[1] == start_coord_alt and i[2] == end_coord:
                middle_ecc = i
            if i[1] == start_coord and i[2] == start_coord_alt:
                middle_ecc = i
    ## if not then reconcile coordinates
        if middle_ecc == 0:
            start_coord_estimate, end_coord_estimate = reconcile_coords(start_coord, end_coord, val)
            for i in val:
                if i[1] == start_coord_estimate and i[2] == end_coord_estimate:
                    middle_ecc = i
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
            if tuple(ecc) in val:
                key_list = list(key)
                key_list[6] = 'yes'
                eccs_perchrom[i] = key_list
    variants_merged = uniq_sort(eccs_perchrom)
    return variants_merged

def ecc_merge(eccs_perchrom):
    variants_grouped = get_variants_grouped(eccs_perchrom)
    representative_variants = get_representative_variants(variants_grouped)
    variants_merged = merge_variants(eccs_perchrom, representative_variants)
    return variants_merged, representative_variants

final_list = []
variants_list = []
## SHOULD NOT HAVE TO DO THIS IN FINAL VERSION OF SCRIPT INPUT SHOULD JUST BE CONFIDENCE_LIST
test_list = [ uniq_sort(confidence_list[i]) for i in range(len(confidence_list))]
for i in range(len(test_list)):
    variants_merged, representative_variants = ecc_merge(test_list[i])
    final_list.append(variants_merged)
    variants_list.append(representative_variants)
flat_list = [item for sublist in final_list for item in sublist]

with open('final.eccs', 'w', newline = '') as final:
    w = csv.writer(final, delimiter = '\t')
    w.writerows(flat_list)

with open('variants_dict', 'w', newline="") as variants_dict:  
    w = csv.writer(variants_dict, delimiter = '\t')
    for i in range(len(variants_list)):
        for key, value in variants_list[i].items():
           w.writerow([key, value])