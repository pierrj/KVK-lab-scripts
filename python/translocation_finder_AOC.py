import csv
import numpy as np
import os
import shutil
from os import listdir
import sys

def get_match_lists(match_file, genomesize_file):
    ref_scaffold_length_dict = {}
    with open(genomesize_file, newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        for row in file_reader:
            ref_scaffold_length_dict[row[0]] = int(row[1])
    match_lists = []
    with open(match_file, newline = '') as file:
            file_reader = csv.reader(file, delimiter = '\t')
            for row in file_reader:
                match_lists.append([row[4], int(row[0]),int(row[1]), row[5], int(row[2]), int(row[3]), ref_scaffold_length_dict[row[4]]])
    ## index and sort
    num_lines = sum(1 for line in open(genomesize_file))
    match_lists_indexed = [[] for i in range(num_lines)]
    scaffold_to_index_dict = {}
    with open(genomesize_file, newline = '') as file:
            file_reader = csv.reader(file, delimiter = '\t')
            for index, scaffold in enumerate(file_reader):
                scaffold_to_index_dict[scaffold[0]] = index
    for i in range(len(match_lists)):
        scaffold_num = scaffold_to_index_dict[match_lists[i][0]]
        match_lists_indexed[scaffold_num].append(match_lists[i])
    match_list_arrays = []
    for i in range(len(match_lists_indexed)):
        match_list_arrays.append(np.array(match_lists_indexed[i], dtype=object))
    return match_list_arrays

def get_candidate_translocations(match_list_processed, query_genomesize_file, tolerance, isclose_percent):
    scaffold_to_index_dict = {}
    with open(query_genomesize_file, newline = '') as file:
            file_reader = csv.reader(file, delimiter = '\t')
            for index, scaffold in enumerate(file_reader):
                scaffold_to_index_dict[scaffold[0]] = index
    scaffold_length_dict = {}
    with open(query_genomesize_file, newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        for row in file_reader:
            scaffold_length_dict[row[0]] = int(row[1])
    translocations = []
    for i in range(len(match_list_processed)):
        matches_per_scaffold = match_list_processed[i]
        for g in range(len(matches_per_scaffold)):
            match = matches_per_scaffold[g]
            start_match = min(match[1], match[2])
            end_match = max(match[1], match[2])
            match_length = end_match - start_match
            match_scaffold = match[3]
            bests = find_translocation_all_possible(start_match, end_match, match_length, matches_per_scaffold, match_scaffold,
                                                scaffold_length_dict, tolerance, isclose_percent)
            if len(bests) > 0:
                for best in bests:
                    best_1 = best[0]
                    best_2 = best[1]
                    translocations.append([match, best_1, best_2])
    return translocations

def find_translocation_all_possible(start_match, end_match, match_length, matches_per_scaffold, match_scaffold, 
                       scaffold_length_dict, tolerance, isclose_percent):
    upstream_matches = matches_per_scaffold[np.logical_and(matches_per_scaffold[:,2] <= start_match+tolerance/2,
                                                          matches_per_scaffold[:,2] >= start_match-tolerance/2)]
    # downstream, starts less than tolerance/2 after end or 20 bp before
    downstream_matches = matches_per_scaffold[np.logical_and(matches_per_scaffold[:,1] >= end_match-tolerance/2,
                                                             matches_per_scaffold[:,1] <= end_match+tolerance/2)]
    bests = []
    for scaffold in scaffold_length_dict.keys():
        if scaffold != match_scaffold:
            scaffold_end = scaffold_length_dict[scaffold]
            subset_upstream = upstream_matches[upstream_matches[:,3] == scaffold]
            subset_downstream = downstream_matches[downstream_matches[:,3] == scaffold]
            # we only care about matches that are longer than the translocation length or start/end at the scaffold borders
            # subset_upstream = subset_upstream[np.logical_or.reduce((subset_upstream[:,4] <= tolerance,
            #                                                subset_upstream[:,5] <= tolerance,
            #                                                subset_upstream[:,4] >= scaffold_end-tolerance,
            #                                                subset_upstream[:,5] >= scaffold_end-tolerance,
            #                                                subset_upstream[:,1] <= tolerance,
            #                                                subset_upstream[:,2] <= tolerance,
            #                                                subset_upstream[:,1] >= subset_upstream[:,6]-tolerance,
            #                                                subset_upstream[:,2] >= subset_upstream[:,6]-tolerance, 
            #                                                abs(subset_upstream[:,4] - subset_upstream[:,5]) >= 1000))]
            # subset_downstream = subset_downstream[np.logical_or.reduce((subset_downstream[:,4] <= tolerance,
            #                                                subset_downstream[:,5] <= tolerance,
            #                                                subset_downstream[:,4] >= scaffold_end-tolerance,
            #                                                subset_downstream[:,5] >= scaffold_end-tolerance,
            #                                                subset_downstream[:,1] <= tolerance,
            #                                                subset_downstream[:,2] <= tolerance,
            #                                                subset_downstream[:,1] >= subset_downstream[:,6]-tolerance,
            #                                                subset_downstream[:,2] >= subset_downstream[:,6]-tolerance, 
            #                                                abs(subset_downstream[:,4] - subset_downstream[:,5]) >= 1000))]
            subset_upstream_plus = subset_upstream[subset_upstream[:,4] - subset_upstream[:,5] < 0]
            subset_downstream_plus = subset_downstream[subset_downstream[:,4] - subset_downstream[:,5] < 0]
            subset_upstream_minus = subset_upstream[subset_upstream[:,4] - subset_upstream[:,5] > 0]
            subset_downstream_minus = subset_downstream[subset_downstream[:,4] - subset_downstream[:,5] > 0]
            for upstream_plus in subset_upstream_plus:
                subset_downstream_plus_subset = subset_downstream_plus[abs(upstream_plus[5]-subset_downstream_plus[:,4]) <= tolerance]
                for downstream_plus in subset_downstream_plus_subset:
                    bests.append([upstream_plus, downstream_plus])
            for upstream_minus in subset_upstream_minus:
                subset_downstream_minus_subset = subset_downstream_minus[abs(upstream_minus[5]-subset_downstream_minus[:,4]) <= tolerance]
                for downstream_minus in subset_downstream_minus_subset:
                    bests.append([upstream_minus, downstream_minus])
    return bests

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

def write_alignments(final_alignments, reference_genomesize, query_genomesize, outdir, reference_name, query_name):
    comparison = reference_name + '_v_' + query_name
    outdir_for_comparison = outdir+'/'+comparison
    os.mkdir(outdir_for_comparison)
    for i in range(len(final_alignments)):
        alignment = final_alignments[i]
        alignment_1 = alignment[0]
        other_chrom_1 = alignment[1]
        other_chrom_2 = alignment[2]
        start_ref = min([alignment_1[1], alignment_1[2]])
        end_ref = max([alignment_1[1], alignment_1[2]])
        start_quer_1 = min([alignment_1[4], alignment_1[5]])
        end_quer_1 = max([alignment_1[4], alignment_1[5]])
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
            w = csv.writer(output_csv, delimiter = '\t') ## lineterminator="\n" might fix issues with windows
            w.writerow([alignment_1[0], start_ref, end_ref])
        with open(outdir_for_comparison + "/" + query_name + '_' + str(i) +".bed", 'w', newline = '') as output_csv:
            w = csv.writer(output_csv, delimiter = '\t') ## lineterminator="\n" might fix issues with windows
            w.writerow([alignment_1[3], start_quer_1, end_quer_1])
            w.writerow([other_chrom_1[3], start_quer_2, end_quer_2])
    return None

processed_matches_file = str(sys.argv[1])
ref_genomesize_file = str(sys.argv[2])
quer_genomesize_file = str(sys.argv[3])
output_directory = str(sys.argv[4])
ref_name = str(sys.argv[5])
quer_name = str(sys.argv[6])

print('started '+quer_name)

match_list = get_match_lists(processed_matches_file, ref_genomesize_file)
translocations = get_candidate_translocations(match_list, quer_genomesize_file, 40,0.1)

if len(translocations) > 0:
    ref_genomesize_dict, quer_genomesize_dict = get_genomesize_dicts(ref_genomesize_file, quer_genomesize_file)
    write_alignments(translocations, ref_genomesize_dict, quer_genomesize_dict, output_directory, ref_name, quer_name)

print('finished '+quer_name)