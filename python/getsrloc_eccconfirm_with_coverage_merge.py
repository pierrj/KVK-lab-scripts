import csv
import ipyparallel as ipp
from itertools import groupby
from itertools import compress
import subprocess
import statistics
import collections


## length filter junctions
with open('samechromosome.exactlytwice.reverseread1.G3_1A_bwamem.bam', newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    with open('getsrloc_test', 'w', newline = '') as filtered:
        w = csv.writer(filtered, delimiter = '\t')
        for line1 in file_reader:
            line2 = next(file_reader)
            sum = int(line1[3]) - int(line2[3])
            if abs(sum) <= 1000000:
                w.writerow(line1)
                w.writerow(line2)

filename = 'getsrloc_test'
filtered_filename = 'qualityfiltered.' + filename
exactlytwice_filename = 'exactlytwice.' + filtered_filename
actualbam_filename = 'actualbam.' + exactlytwice_filename
bedfile = filename + '.bed'
samtools = '/global/home/groups/consultsw/sl-7.x86_64/modules/samtools/1.8/bin/samtools'
bedtools = '/global/home/groups/consultsw/sl-7.x86_64/modules/bedtools/2.28.0/bin/bedtools'

quality_filter = '''awk '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\\\1", "", $6); if((a !~ /[DMIHS]/ && int(a) > 49 ) || (b !~ /[DMIHS]/ && int(b) > 49)) print $0}' ''' + filename + ''' > ''' + filtered_filename
exactlytwice_filter = '''awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' ''' + filtered_filename + ''' ''' + filtered_filename + ''' > ''' + exactlytwice_filename
make_bam = 'bash -c \"' + samtools + ' view -b -h <(cat <(' + samtools + ' view -H G3_1A_bwamem.bam) ' + exactlytwice_filename +') > ' + actualbam_filename + '\"'
bamtobed_sort = bedtools + ' bamtobed -i ' + actualbam_filename + ' | sort -k 4,4 -k 2,2 > ' + bedfile

subprocess.run(quality_filter, shell= True)
subprocess.run(exactlytwice_filter, shell= True)
subprocess.run(make_bam, shell= True)
subprocess.run(bamtobed_sort, shell= True)

with open('getsrloc_test.bed', newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    with open('merged.getsrloc_test', 'w', newline = '') as filtered:
        w = csv.writer(filtered, delimiter = '\t')
        for line1 in file_reader:
            line2 = next(file_reader)
            line1[2] = line2[1]
            w.writerow(line1)

with open('merged.getsrloc_test', newline = '') as eccloc:
    eccloc_reader = csv.reader(eccloc, delimiter = '\t')
    eccloc_list = [[int(row[0][10:12]) - 1, int(row[1]), int(row[2])] for row in eccloc_reader]
    
with open('discordantmappedreads.oppositefacing.bed', newline = '') as discordant:
    discordant_reader = csv.reader(discordant, delimiter = '\t')
    discordant_indexed = [[] for i in range(56)]
    for row in discordant_reader:
        discordant_indexed[(int(row[0][10:12])-1)].append([int(row[1]), int(row[2])])

def confirmeccs(ecc):
    for i in range(0, len(discordant_indexed[ecc[0]]), 2):
        reada = discordant_indexed[ecc[0]][i]
        readb = discordant_indexed[ecc[0]][i+1]
        if reada[0] < readb[0]:
            read1 = reada
            read2 = readb
        else:
            read1 = readb
            read2 = reada
        if ecc[1] <= read1[0] <= ecc[2] and ecc[1] <= read1[1] <= ecc[2] and ecc[1] <= read2[0] <= ecc[2] and ecc[1] <= read2[1] <= ecc[2]:
            distance = abs(read1[1] - ecc[1]) + abs(read2[0] - ecc[2])
            if distance <= 500:
                return True
    return False

rc = ipp.Client(profile='default', cluster_id='')
dview = rc[:]
dview.block = True
lview = rc.load_balanced_view()
lview.block = True

mydict = dict(discordant_indexed = discordant_indexed)
dview.push(mydict)

yesornoeccs = list(lview.map(confirmeccs, eccloc_list))
confirmedeccs = list(compress(eccloc_list, yesornoeccs))

with open('parallel.confirmed', 'w', newline = '') as confirmed:
    w = csv.writer(confirmed, delimiter = '\t')
    w.writerows(confirmedeccs)

with open('parallel.confirmed') as confirmed:
    confirmed_reader = csv.reader(confirmed, delimiter = '\t')
    confirmed_list = [[int(row[0]), int(row[1]), int(row[2])] for row in confirmed_reader]

def splitreadcount_filter(lst, k):
    tuple_list = [tuple(x) for x in lst]
    counted = collections.Counter(tuple_list)
    count_filtered = [el for el in lst if counted[tuple(el)] >= k]
    return [list(x) for x in set(tuple(x) for x in count_filtered)]

confirmed_list = splitreadcount_filter(confirmed_list, 1)

with open('genomecoverage.mergedandpe.G3_1A_bwamem.bed') as coverage:
    coverage_reader = csv.reader(coverage, delimiter = '\t')
    coverage_indexed = [[] for i in range(56)]
    for row in coverage_reader:
        coverage_indexed[(int(row[0][10:12])-1)].append([int(row[1]) -1, int(row[2])])

def confidence_check(confirmed):
    for i in range(len(confirmed)):
        ecc = confirmed[i]
        region_len = ecc[2] - ecc[1]
        region = coverage_indexed[ecc[0]][ecc[1]:(ecc[2]+1)]
        region_cov = [region[k][1] for k in range(len(region))]
        mean_region = round(statistics.mean(region_cov), 2)
        beforestart = ecc[1] - 1 - region_len
        afterstart = ecc[2] + 1
        region_before = coverage_indexed[ecc[0]][beforestart:beforestart + 1 + region_len]
        region_before_cov = [region_before[j][1] for j in range(len(region_before))]
        region_after = coverage_indexed[ecc[0]][afterstart:afterstart + 1 + region_len]
        region_after_cov = [region_after[g][1] for g in range(len(region_after))]
        if beforestart > 0:
            mean_before = round(statistics.mean(region_before_cov), 2)
        else:
            mean_before = 'N/A'
        if afterstart + region_len <= coverage_indexed[ecc[0]][-1][0]:
            mean_after = round(statistics.mean(region_after_cov), 2)
        else:
            mean_after = 'N/A'
        coverage_string = str(mean_before)+';'+str(mean_region)+';'+str(mean_after)
        if region_cov.count(0) / len(region_cov) > 0.05:
            ecc.append('lowq')
            ecc.append('incomplete_coverage')
            ecc.append(coverage_string)
            continue
        if beforestart <= 0:
            ecc.append('conf')
            ecc.append('too_close_to_start')
            ecc.append(coverage_string)
            continue
        if afterstart + region_len > coverage_indexed[ecc[0]][-1][0]:
            ecc.append('conf')
            ecc.append('too_close_to_end')
            ecc.append(coverage_string)
            continue
        if mean_region >= (2*mean_before) and mean_region >= (2*mean_after):
            ecc.append('hconf')
            ecc.append('none')
            ecc.append(coverage_string)
        else:
            ecc.append('conf')
            ecc.append('coverage_too_low')
            ecc.append(coverage_string)
    return confirmed

confidence_list = confidence_check(confirmed_list)

def index_confidence_list(confidence_list):
    confidence_index = [[] for i in range(56)] ## scaffold number here, get this variable from somewhere else
    for i in range(len(confidence_list)):
        row = confidence_list[i]
        confidence_index[int(row[0])].append([int(row[0])+1, int(row[1]), int(row[2]), str(row[3]), str(row[4]), str(row[5]), 'no']) ## adds no to be changed to yes if their is a variant
    return confidence_index

confidence_indexed = index_confidence_list(confidence_list)

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
#test_list = [ sorted(confidence_indexed[i], key=lambda x:x [1]) for i in range(len(confidence_indexed))]
for i in range(len(confidence_indexed)):
    variants_merged, representative_variants = ecc_merge(confidence_indexed[i])
    final_list.append(variants_merged)
    variants_list.append(representative_variants)
flat_list = [item for sublist in final_list for item in sublist]

with open('ecccaller_output.details.tsv', 'w', newline = '') as final:
    w = csv.writer(final, delimiter = '\t')
    w.writerows(flat_list)

with open('ecccaller_output.bed', 'w', newline = '') as bed:
    w = csv.writer(bed, delimiter = '\t')
    for i in range(len(flat_list)):
        scaffold_string = 'MQOP010000' + str(flat_list[i][0]).zfill(2) + ".1"
        row = [scaffold_string, flat_list[i][1], flat_list[i][2]]
        w.writerow(row)

with open('ecccaller_variants.tsv', 'w', newline="") as variants_dict:
    w = csv.writer(variants_dict, delimiter = '\t')
    for i in range(len(variants_list)):
        for key, value in variants_list[i].items():
           w.writerow([key, value])