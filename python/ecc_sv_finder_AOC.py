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
import os
import shutil
from os import listdir
import sys

## USAGE ##
# this script parses nucmer alignments and look for a specific set of alignments that are indicative of ecc mediated translocation
# "processed_matches_file" coords file as input
# "ref_genomesize_file" file containing scaffold names and scaffold lengths for reference
# "quer_genomesize_file" file containing scaffold names and scaffold lengths for querry
# "output_directory" directory for outputs
# "ref_name" reference prefix/name
# "quer_name" query prefix/name

processed_matches_file = str(sys.argv[1])
ref_genomesize_file = str(sys.argv[2])
quer_genomesize_file = str(sys.argv[3])
output_directory = str(sys.argv[4])
ref_name = str(sys.argv[5])
quer_name = str(sys.argv[6])

def get_match_lists(match_file, genomesize_file):
    ref_scaffold_length_dict = {} # read in genomesize file
    with open(genomesize_file, newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        for row in file_reader:
            ref_scaffold_length_dict[row[0]] = int(row[1])
    match_lists = [] # read in coords file
    with open(match_file, newline = '') as file:
            file_reader = csv.reader(file, delimiter = '\t')
            for row in file_reader:
                match_lists.append([row[4], int(row[0]),int(row[1]), row[5], int(row[2]), int(row[3]), ref_scaffold_length_dict[row[4]]])
    ## index and sort
    num_lines = sum(1 for line in open(genomesize_file))
    match_lists_indexed = [[] for i in range(num_lines)]
    scaffold_to_index_dict = {} # this is used to order scaffolds and retrieve them by their older later
    with open(genomesize_file, newline = '') as file:
            file_reader = csv.reader(file, delimiter = '\t')
            for index, scaffold in enumerate(file_reader):
                scaffold_to_index_dict[scaffold[0]] = index
    for i in range(len(match_lists)):
        scaffold_num = scaffold_to_index_dict[match_lists[i][0]]
        match_lists_indexed[scaffold_num].append(match_lists[i])
    match_list_arrays = [] # turn match list into arrays for numpy
    for i in range(len(match_lists_indexed)):
        match_list_arrays.append(np.array(match_lists_indexed[i], dtype=object))
    return match_list_arrays

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
    overlap_percent = np.vstack((overlap_or_zero/lst_length, overlap_or_zero/twod_array_length)).T.max(axis=1)
    boolean_array = overlap_percent < percent # compare to percentage cutoff
    return boolean_array

def get_promising_alignments(match_list_processed,tolerance): # get a pair of alignments that are nearby in query, which indicates a translocation event
    promising_alignments = []
    for i in range(len(match_list_processed)):
        scaffold = match_list_processed[i]
        for g in range(len(scaffold)): # same scaffold
            match = scaffold[g]
            ## find another match on the same reference scaffold that is nearby
            ## either end is within tolerance of end
            ## or start is within tolerance of end
            nearby_matches = scaffold[np.logical_or(
                        abs(scaffold[:,1] - match[2]) <= tolerance,
                        abs(scaffold[:,2] - match[1]) <= tolerance)]
            match_same_chrom = nearby_matches[nearby_matches[:,3] == match[3]]
            nearby_on_same_chrom = match_same_chrom[np.logical_or(abs(match_same_chrom[:,4] - match[5]) <= tolerance,
                                                                   abs(match_same_chrom[:,5] - match[4]) <= tolerance)]
            if nearby_on_same_chrom.size > 0: # check if any are left
                nearby_on_same_chrom = nearby_on_same_chrom[np.logical_and(nearby_on_same_chrom[:,1] != match[1],
                                                                          nearby_on_same_chrom[:,1] != match[2])] # no overlaps
                ## this is written a little awkwardly, probably could've been if statements instead
                ## but basically, one of the two is always empty here
                ## so basically it is one or the other, i just wrote both instead of an if statement
                nearby_on_same_chrom_match_is_first = nearby_on_same_chrom[nearby_on_same_chrom[:,1:2].min(axis=1) > min([match[1], match[2]])] # get relative location
                nearby_on_same_chrom_match_is_second = nearby_on_same_chrom[nearby_on_same_chrom[:,1:2].min(axis=1) < min([match[1], match[2]])]
                if match[4] - match[5] > 0:
                    match_order = '-'
                elif match[4] - match[5] < 0 :
                    match_order = '+'
                if match_order == '+': ## now looking at query matches corresponding to identified matches to see if they are in a proper order and orientation
                    nearby_on_same_chrom_match_is_first = nearby_on_same_chrom_match_is_first[nearby_on_same_chrom_match_is_first[:,4] - nearby_on_same_chrom_match_is_first[:,5] < 0]
                    nearby_on_same_chrom_match_is_second = nearby_on_same_chrom_match_is_second[nearby_on_same_chrom_match_is_second[:,4] - nearby_on_same_chrom_match_is_second[:,5] < 0]
                    nearby_on_same_chrom_match_is_first = nearby_on_same_chrom_match_is_first[nearby_on_same_chrom_match_is_first[:,4:6].max(axis=1) < max([match[4], match[5]])]
                    nearby_on_same_chrom_match_is_second = nearby_on_same_chrom_match_is_second[nearby_on_same_chrom_match_is_second[:,4:6].max(axis=1) > max([match[4], match[5]])]
                elif match_order == '-':
                    nearby_on_same_chrom_match_is_first = nearby_on_same_chrom_match_is_first[nearby_on_same_chrom_match_is_first[:,4] - nearby_on_same_chrom_match_is_first[:,5] > 0]
                    nearby_on_same_chrom_match_is_second = nearby_on_same_chrom_match_is_second[nearby_on_same_chrom_match_is_second[:,4] - nearby_on_same_chrom_match_is_second[:,5] > 0]
                    nearby_on_same_chrom_match_is_first = nearby_on_same_chrom_match_is_first[nearby_on_same_chrom_match_is_first[:,4:6].max(axis=1) > max([match[4], match[5]])]
                    nearby_on_same_chrom_match_is_second = nearby_on_same_chrom_match_is_second[nearby_on_same_chrom_match_is_second[:,4:6].max(axis=1) < max([match[4], match[5]])]
                # making sure the alignments dont overlap in query or reference by more than 20%
                nearby_on_same_chrom_match_is_first = nearby_on_same_chrom_match_is_first[numpy_get_overlap_percent([match[1], match[2]], nearby_on_same_chrom_match_is_first[:,1:3] ,0.2)]
                nearby_on_same_chrom_match_is_first = nearby_on_same_chrom_match_is_first[numpy_get_overlap_percent([match[4], match[5]], nearby_on_same_chrom_match_is_first[:,4:6] ,0.2)]
                nearby_on_same_chrom_match_is_second = nearby_on_same_chrom_match_is_second[numpy_get_overlap_percent([match[1], match[2]], nearby_on_same_chrom_match_is_second[:,1:3] ,0.2)]
                nearby_on_same_chrom_match_is_second = nearby_on_same_chrom_match_is_second[numpy_get_overlap_percent([match[4], match[5]], nearby_on_same_chrom_match_is_second[:,4:6] ,0.2)]
                # write promising matches as lists
                ## so far we have only found pairs of alignments that are close together on the reference scaffold and properly oriented on the query scaffold
                ## this likely represents some sort of translocation, as long as there is an alignment in the middle that is from a different scaffold
                for ecc_match in nearby_on_same_chrom_match_is_first:
                    if [match.tolist(), ecc_match.tolist()] not in promising_alignments:
                        promising_alignments.append([match.tolist(), ecc_match.tolist()])
                for ecc_match in nearby_on_same_chrom_match_is_second:
                    if [ecc_match.tolist(), match.tolist()] not in promising_alignments:
                        promising_alignments.append([ecc_match.tolist(), match.tolist()])
    return promising_alignments

def get_final_alignments(promising_alignments, match_list, ref_genomesize_file, query_genomesize_file, tolerance, isclose_percent): # process promising alignments, looking for alignments that are found between the pairs
    # process genome size files
    scaffold_to_index_dict = {}
    with open(ref_genomesize_file, newline = '') as file:
            file_reader = csv.reader(file, delimiter = '\t')
            for index, scaffold in enumerate(file_reader):
                scaffold_to_index_dict[scaffold[0]] = index
    scaffold_length_dict = {}
    with open(query_genomesize_file, newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        for row in file_reader:
            scaffold_length_dict[row[0]] = int(row[1])
    # get translocations
    translocations = []
    for i in range(len(promising_alignments)):
        alignment = promising_alignments[i]
        alignment_1 = alignment[0]
        alignment_2 = alignment[1]
        start_match = min([alignment_1[1], alignment_1[2], alignment_2[1], alignment_2[2]])
        end_match = max([alignment_1[1], alignment_1[2], alignment_2[1], alignment_2[2]])
        match_length = end_match - start_match
        matches_per_scaffold = match_list[scaffold_to_index_dict[alignment_1[0]]] # matches per scaffold
        match_scaffold = alignment_1[3]
        ## get translocations using this function
        bests = find_translocation_all_possible(start_match, end_match, match_length, matches_per_scaffold, match_scaffold,
                                            scaffold_length_dict, tolerance, isclose_percent)
        # write all possible translocations
        # at this point we should have something that looks like an ecc sv
        # this means, a translocation 
        ### two alignments, next to each other in the reference but with a gap in between them in the query
        # but with a sign of eccdna mediated sv
        ### the gap in between, rather than being a continuous alignment, is actually two alignments and their order is swapped (see figure examples for a visualization)
        if len(bests) > 0:
            for best in bests:
                best_1 = best[0]
                best_2 = best[1]
                translocations.append([alignment_1, alignment_2, best_1, best_2])
    return translocations

def find_translocation_all_possible(start_match, end_match, match_length, matches_per_scaffold, match_scaffold, 
                       scaffold_length_dict, tolerance, isclose_percent):
    # get matches downstream of end of alignment 1
    # and get matches upstream of start of alignment 2
    ## basically to get alignments sandwhiched between promising alignments
    upstream_matches = matches_per_scaffold[np.logical_and(matches_per_scaffold[:,2] <= start_match+tolerance/2,
                                                          matches_per_scaffold[:,2] >= start_match-tolerance/2)]
    downstream_matches = matches_per_scaffold[np.logical_and(matches_per_scaffold[:,1] >= end_match-tolerance/2,
                                                             matches_per_scaffold[:,1] <= end_match+tolerance/2)]
    bests = []
    for scaffold in scaffold_length_dict.keys():
        if scaffold != match_scaffold: # make sure alignments are not the same scaffold as the original pair of matches
            # look through matches by scaffold
            subset_upstream = upstream_matches[upstream_matches[:,3] == scaffold]
            subset_downstream = downstream_matches[downstream_matches[:,3] == scaffold]
            # to find pairs of matches that are properly oriented for an ecc sv in the query
            ## again this is a reversed order from what you would expect from a usual translocation
            subset_upstream_plus = subset_upstream[subset_upstream[:,4] - subset_upstream[:,5] < 0]
            subset_downstream_plus = subset_downstream[subset_downstream[:,4] - subset_downstream[:,5] < 0]
            subset_upstream_minus = subset_upstream[subset_upstream[:,4] - subset_upstream[:,5] > 0]
            subset_downstream_minus = subset_downstream[subset_downstream[:,4] - subset_downstream[:,5] > 0]
            # make sure they are close enough and write
            for upstream_plus in subset_upstream_plus:
                subset_downstream_plus_subset = subset_downstream_plus[abs(upstream_plus[5]-subset_downstream_plus[:,4]) <= tolerance]
                for downstream_plus in subset_downstream_plus_subset:
                    bests.append([upstream_plus, downstream_plus])
            for upstream_minus in subset_upstream_minus:
                subset_downstream_minus_subset = subset_downstream_minus[abs(upstream_minus[5]-subset_downstream_minus[:,4]) <= tolerance]
                for downstream_minus in subset_downstream_minus_subset:
                    bests.append([upstream_minus, downstream_minus])
    return bests

# read in genome size files...
def get_genomesize_dicts(reference_genomesize, query_genomesize):
    reference_genomesize_dict = {}
    query_genomesize_dict = {}
    with open(reference_genomesize, newline = '') as file:
            file_reader = csv.reader(file, delimiter = '\t')
            for row in file_reader:
                reference_genomesize_dict[row[0]] = int(row[1])
    with open(query_genomesize, newline = '') as file:
            file_reader = csv.reader(file, delimiter = '\t')
            for row in file_reader:
                query_genomesize_dict[row[0]] = int(row[1])
    return reference_genomesize_dict, query_genomesize_dict

## finally write fasta files to align to each other, making sure not to go past scaffold edges
def write_alignments(final_alignments, reference_genomesize, query_genomesize, outdir, reference_name, query_name):
    comparison = reference_name + '_v_' + query_name
    outdir_for_comparison = outdir+'/'+comparison
    os.mkdir(outdir_for_comparison)
    for i in range(len(final_alignments)):
        alignment = final_alignments[i]
        alignment_1 = alignment[0]
        alignment_2 = alignment[1]
        other_chrom_1 = alignment[2]
        other_chrom_2 = alignment[3]
        start_ref = min([alignment_1[1], alignment_1[2], alignment_2[1], alignment_2[2]])
        end_ref = max([alignment_1[1], alignment_1[2], alignment_2[1], alignment_2[2]])
        start_quer_1 = min([alignment_1[4], alignment_1[5], alignment_2[4], alignment_2[5]])
        end_quer_1 = max([alignment_1[4], alignment_1[5], alignment_2[4], alignment_2[5]])
        ## just making sure to get a sequence with 20k bp on either side
        ## but also making sure not to go past either side of the scaffolds
        if other_chrom_1[4] - other_chrom_1[5] > 0:
            match_order = '-'
        elif other_chrom_1[4] - other_chrom_1[5] < 0 :
            match_order = '+'
        if match_order == '-':
            start_quer_2 = other_chrom_2[4]
            end_quer_2 = other_chrom_1[5]
        if match_order == '+':
            start_quer_2 = other_chrom_1[5]
            end_quer_2 = other_chrom_2[4]
        if start_ref - 20000 > 0:
            start_ref = start_ref - 20000
        else:
            start_ref = 1
        if end_ref + 20000 < reference_genomesize[alignment_1[0]]:
            end_ref = end_ref + 20000
        else:
            end_ref = reference_genomesize[alignment_1[0]]
        if start_quer_1 - 20000 > 0:
            start_quer_1 = start_quer_1 - 20000
        else:
            start_quer_1 = 1
        if end_quer_1 + 20000 < query_genomesize[alignment_1[3]]:
            end_quer_1 = end_quer_1 + 20000
        else:
            end_quer_1 = query_genomesize[alignment_1[3]]
        if start_quer_2 - 20000 > 0:
            start_quer_2 = start_quer_2 - 20000
        else:
            start_quer_2 = 1
        if end_quer_2 + 20000 < query_genomesize[other_chrom_1[3]]:
            end_quer_2 = end_quer_2 + 20000
        else:
            end_quer_2 = query_genomesize[other_chrom_1[3]]
        with open(outdir_for_comparison + "/" + reference_name + '_' +str(i) +".bed", 'w', newline = '') as output_csv:
            w = csv.writer(output_csv, delimiter = '\t')
            w.writerow([alignment_1[0], start_ref, end_ref])
        with open(outdir_for_comparison + "/" + query_name + '_' + str(i) +".bed", 'w', newline = '') as output_csv:
            w = csv.writer(output_csv, delimiter = '\t')
            w.writerow([alignment_1[3], start_quer_1, end_quer_1])
            w.writerow([other_chrom_1[3], start_quer_2, end_quer_2])
    return None

print('started '+quer_name)

# run everything
match_list = get_match_lists(processed_matches_file, ref_genomesize_file)
promising_list = get_promising_alignments(match_list, 100)
final_list = get_final_alignments(promising_list, match_list, ref_genomesize_file, quer_genomesize_file, 100, 0.1)

# only write files/make directories if any were found
if len(final_list) > 0:
    ref_genomesize_dict, quer_genomesize_dict = get_genomesize_dicts(ref_genomesize_file, quer_genomesize_file)
    write_alignments(final_list, ref_genomesize_dict, quer_genomesize_dict, output_directory, ref_name, quer_name)

print('finished '+quer_name)