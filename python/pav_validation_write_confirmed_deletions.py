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
import sys
import more_itertools as mit
from itertools import filterfalse
from os import listdir

## USAGE ##
# this script parses alignments of putative gene deletion locations in guy11 vs other genomes and makes sure that the gap in the alignment overlaps with the location of the gene
# "directory_coords" directory of nucmer outputs
# "directory_gene_locs" directory of gene locations in the alignment
# "output" output files, list of gene deletion locations in the guy11 genome

directory_coords = str(sys.argv[1])
directory_gene_locs = str(sys.argv[2])
output = str(sys.argv[3])

## yields runs of consecutive numbers
def find_ranges(iterable):
    for group in mit.consecutive_groups(iterable):
        group = list(group)
        if len(group) == 1:
            yield group[0], group[0]
        else:
            yield group[0], group[-1]

# needs to be sorted to work on linux...
coords_directory_list = sorted(listdir(directory_coords))
genelocs_directory_list = sorted(listdir(directory_gene_locs))

confirmed_deletions = []
unconfirmed_deletions = []
for i in range(len(coords_directory_list)): ## loop through coords directory
    coord = coords_directory_list[i]
    gene_loc = genelocs_directory_list[i]
    gene = gene_loc[:13] # get gene location in alignment
    genome = coord[14:-24]
    with open(directory_coords+'/'+coord, newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        deletion_processed_coords_list = []
        for row in file_reader:
            deletion_processed_coords_list.append([row[4], int(row[0]), int(row[1]), row[5], int(row[2]), int(row[3])]) ## read in coords files
    with open(directory_gene_locs+'/'+gene_loc, newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        for row in file_reader:
            gene_loc_list = [row[0], int(row[1]), int(row[2])] # gene locations
    start_region = int(deletion_processed_coords_list[0][0].split(':')[1].split('-')[0]) # get alignment regions
    end_region = int(deletion_processed_coords_list[0][0].split(':')[1].split('-')[1])
    region_length = end_region-start_region
    current_range = range(1,region_length) ## get regions where alignments are missing based off ranges where alignments are presents
    for match in deletion_processed_coords_list:
        current_range = list(filterfalse(lambda x: match[1] <= x <= match[2], current_range))
    deletions = list(find_ranges(current_range))
    matching_deletion = []
    start_gene = gene_loc_list[1]
    end_gene = gene_loc_list[2]
    for deletion in deletions: ## see if deletion overlaps with gene to any extent
        start_deletion = deletion[0]
        end_deletion = deletion[1]
        if (start_deletion > start_gene and start_deletion < end_gene ) or (end_deletion > start_gene and end_deletion < end_gene) or (start_deletion < start_gene and end_deletion > end_gene):
            matching_deletion.append(deletion)
    for deletion in matching_deletion:
        confirmed_deletions.append([gene_loc_list[0], deletion[0]+start_region, deletion[1]+start_region, gene, genome])

# write confirmed deletions
with open(output , 'w', newline = '') as output_csv:
    w = csv.writer(output_csv, delimiter = '\t')
    for deletion in confirmed_deletions:
        w.writerow([deletion[0], deletion[1], deletion[2], deletion[3], deletion[4]])