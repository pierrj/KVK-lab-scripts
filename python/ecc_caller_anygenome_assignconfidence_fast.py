import csv
import statistics
import collections

coverage_file = str(sys.argv[1])

output_name = str(sys.argv[2])

scaffold_number = int(sys.argv[3])

with open('parallel.confirmed') as confirmed:
    confirmed_reader = csv.reader(confirmed, delimiter = '\t')
    confirmed_list = [[int(row[0]), int(row[1]), int(row[2])] for row in confirmed_reader]

def add_splitreadcount(lst):
    tuple_list = [tuple(x) for x in lst]
    counted = collections.Counter(tuple_list)
    count_filtered = [[el[0], el[1], el[2], counted[tuple(el)]] for el in lst]
    return [list(x) for x in set(tuple(x) for x in count_filtered)]

confirmed_list_sr = add_splitreadcount(confirmed_list)

def index_ecc_list(ecc_list):
    ecc_index = [[] for i in range(scaffold_number)]
    for i in range(len(ecc_list)):
        row = ecc_list[i]
        # add no to be changed to yes later by ecc_merge if variants are found
        ecc_index[int(row[0])].append([int(row[0]), int(row[1]), int(row[2]), int(row[3]), 'no'])
    return ecc_index

ecc_indexed = index_ecc_list(confirmed_list_sr)

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

# get ideal representative eccDNA based off start and end medians of variants
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

# get closest matching coordinate based off one set/forced coordinate
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

# test which coordinate pair is closest to the ideal coordinates
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

# get single representative variants for eccs with multiple variants
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

# get only unique eccDNAs and sort based off start coordinate
def uniq_sort(start_list):
    uniq_list = [list(x) for x in set(tuple(x) for x in start_list)]
    end_list = sorted(uniq_list,key=lambda x: x[1])
    return end_list

# loop through eccDNA list, if the eccDNA is found in any of the values in the representative_variants dictionary
def merge_variants(eccs_perchrom, representative_variants):
    for i in range(len(eccs_perchrom)):
        ecc = eccs_perchrom[i]
        for key, val in representative_variants.items():
            if tuple(ecc) in val and [tuple(ecc)] != val: #double check to make sure val isnt just ecc
                key_list = list(key)
                # replace that ecc with the representative variant key and change the no in the variant column to yes
                key_list[4] = 'yes'
                # replace the number of split reads by the sum of split reads of the eccDNAs
                sr_sum = 0
                for l in range(len(val)):
                    sr_sum += val[l][3]
                key_list[3] = sr_sum
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

merged_list = []
splitreads_list = []
for i in range(len(ecc_indexed)):
    variants_merged, representative_variants = ecc_merge(ecc_indexed[i])
    merged_list.append(variants_merged)
    splitreads_list.append(representative_variants)

flat_merged_list = [item for sublist in merged_list for item in sublist]

with open(coverage_file) as coverage:
    coverage_reader = csv.reader(coverage, delimiter = '\t')
    # index coverage file so that confidence check only looks at the same chromosome
    coverage_indexed = [[] for i in range(scaffold_number)]
    for row in coverage_reader:
        coverage_indexed[(int(row[0])-1)].append([int(row[1]), int(row[2])])

def confidence_check(confirmed):
    for i in range(len(confirmed)):
        ecc = confirmed[i]
        # get coverage of confirmed ecc region
        region_len = ecc[2] - ecc[1]
        region = coverage_indexed[ecc[0]][ecc[1]:(ecc[2]+1)]
        region_cov = [region[k][1] for k in range(len(region))]
        mean_region = round(statistics.mean(region_cov), 2)
        # define regions before and after confirmed ecc that are of the same length as the ecc
        beforestart = ecc[1] - 1 - region_len
        afterstart = ecc[2] + 1
        region_before = coverage_indexed[ecc[0]][beforestart:beforestart + 1 + region_len]
        region_before_cov = [region_before[j][1] for j in range(len(region_before))]
        region_after = coverage_indexed[ecc[0]][afterstart:afterstart + 1 + region_len]
        region_after_cov = [region_after[g][1] for g in range(len(region_after))]
        # check to see if the before and after regions go beyond the length of the chromosome if so dont calculate mean coverage for that region
        if beforestart > 0:
            mean_before = round(statistics.mean(region_before_cov), 2)
        else:
            mean_before = 'N/A'
        if afterstart + region_len <= coverage_indexed[ecc[0]][-1][0]:
            mean_after = round(statistics.mean(region_after_cov), 2)
        else:
            mean_after = 'N/A'
        # write coverage string containing the means of the regions before, within, and after the confirmed ecc
        coverage_string = str(mean_before)+';'+str(mean_region)+';'+str(mean_after)
        # if less than 95% of the ecc region is covered than the ecc is low confidence
        if region_cov.count(0) / len(region_cov) > 0.05:
            ecc.append('lowq')
            ecc.append('incomplete_coverage')
            ecc.append(coverage_string)
            continue
        if ecc[3] == 1:
            ecc.append('lowq')
            ecc.append('only_one_splitread')
            ecc.append(coverage_string)
            continue
        # if the ecc before or after regions fall beyond the borders of the chromosome than the ecc is medium confidence
        if beforestart <= 0 and ecc[3] < 3:
            ecc.append('conf')
            ecc.append('too_close_to_start')
            ecc.append(coverage_string)
            continue
        if afterstart + region_len > coverage_indexed[ecc[0]][-1][0] and ecc[3] < 3:
            ecc.append('conf')
            ecc.append('too_close_to_end')
            ecc.append(coverage_string)
            continue
        if beforestart <= 0 and ecc[3] >= 3:
            ecc.append('hconf')
            ecc.append('splitreads')
            ecc.append(coverage_string)
            continue
        if afterstart + region_len > coverage_indexed[ecc[0]][-1][0] and ecc[3] >= 3:
            ecc.append('hconf')
            ecc.append('splitreads')
            ecc.append(coverage_string)
            continue
        # if the mean coverage of the eccDNA is twice the size of the before and after regions OR the eccDNA is supported by three or more splitreads then eccDNA is high confidence, otherwise is it medium confidence
        if mean_region >= (2*mean_before) and mean_region >= (2*mean_after) and ecc[3] < 3:
            ecc.append('hconf')
            ecc.append('coverage')
            ecc.append(coverage_string)
            continue
        if (mean_region < (2*mean_before) or mean_region < (2*mean_after)) and ecc[3] >= 3:
            ecc.append('hconf')
            ecc.append('splitreads')
            ecc.append(coverage_string)
            continue
        if mean_region >= (2*mean_before) and mean_region >= (2*mean_after) and ecc[3] >= 3:
            ecc.append('hconf')
            ecc.append('splitreads_and_coverage')
            ecc.append(coverage_string)
        else:
            ecc.append('conf')
            ecc.append('coverage_and_splitreads_too_low')
            ecc.append(coverage_string)
    return confirmed

confidence_flat_merged_list = confidence_check(flat_merged_list)

with open('ecccaller_output.' + output_name + '.details.tsv', 'w', newline = '') as bed:
    w = csv.writer(bed, delimiter = '\t')
    for i in range(len(confidence_flat_merged_list)):
        row = [confidence_flat_merged_list[i][0]+1, confidence_flat_merged_list[i][1], confidence_flat_merged_list[i][2], confidence_flat_merged_list[i][3], confidence_flat_merged_list[i][4],confidence_flat_merged_list[i][5], confidence_flat_merged_list[i][6], confidence_flat_merged_list[i][7]]
        w.writerow(row)

# write file with split reads per ecc, with representative ecc in the first column and all split reads in second column
with open('ecccaller_splitreads.' + output_name + '.tsv', 'w', newline="") as variants_dict:
    w = csv.writer(variants_dict, delimiter = '\t')
    for i in range(len(splitreads_list)):
        for key, value in splitreads_list[i].items():
           w.writerow([key, value])

# write bed file output, with actual chromosome names, start and end coordinates and ecc number
with open('ecccaller_output.' + output_name + '.bed', 'w', newline = '') as bed:
    w = csv.writer(bed, delimiter = '\t')
    for i in range(len(confidence_flat_merged_list)):
        if confidence_flat_merged_list[i][5] == 'lowq':
            color = '255,0,0'
        if confidence_flat_merged_list[i][5] == 'conf':
            color = '255,255,0'
        if confidence_flat_merged_list[i][5] == 'hconf':
            color = '0,255,0'
        row = [confidence_flat_merged_list[i][0]+1, confidence_flat_merged_list[i][1], confidence_flat_merged_list[i][2], i, 0, '+', confidence_flat_merged_list[i][1], confidence_flat_merged_list[i][2], color ]
        w.writerow(row)