import csv
import sys
import more_itertools as mit
from itertools import filterfalse
from os import listdir

directory_coords = str(sys.argv[1])
directory_gene_locs = str(sys.argv[2])
output = str(sys.argv[3])

def find_ranges(iterable):
    """Yield range of consecutive numbers."""
    for group in mit.consecutive_groups(iterable):
        group = list(group)
        if len(group) == 1:
            yield group[0], group[0]
        else:
            yield group[0], group[-1]

coords_directory_list = sorted(listdir(directory_coords))
genelocs_directory_list = sorted(listdir(directory_gene_locs))

confirmed_deletions = []
unconfirmed_deletions = []
for i in range(len(coords_directory_list)):
    coord = coords_directory_list[i]
    gene_loc = genelocs_directory_list[i]
    gene = gene_loc[:13]
    genome = coord[14:-24]
    with open(directory_coords+'/'+coord, newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        deletion_processed_coords_list = []
        for row in file_reader:
            deletion_processed_coords_list.append([row[4], int(row[0]), int(row[1]), row[5], int(row[2]), int(row[3])])
    with open(directory_gene_locs+'/'+gene_loc, newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        for row in file_reader:
            gene_loc_list = [row[0], int(row[1]), int(row[2])]
    start_region = int(deletion_processed_coords_list[0][0].split(':')[1].split('-')[0])
    end_region = int(deletion_processed_coords_list[0][0].split(':')[1].split('-')[1])
    region_length = end_region-start_region
    current_range = range(1,region_length)
    for match in deletion_processed_coords_list:
        current_range = list(filterfalse(lambda x: match[1] <= x <= match[2], current_range))
    deletions = list(find_ranges(current_range))
    matching_deletion = []
    start_gene = gene_loc_list[1]
    end_gene = gene_loc_list[2]
    for deletion in deletions:
        start_deletion = deletion[0]
        end_deletion = deletion[1]
        if (start_deletion > start_gene and start_deletion < end_gene ) or (end_deletion > start_gene and end_deletion < end_gene) or (start_deletion < start_gene and end_deletion > end_gene):
            matching_deletion.append(deletion)
    for deletion in matching_deletion:
        confirmed_deletions.append([gene_loc_list[0], deletion[0]+start_region, deletion[1]+start_region, gene, genome])

with open(output , 'w', newline = '') as output_csv:
    w = csv.writer(output_csv, delimiter = '\t') ## lineterminator="\n" might fix issues with windows
    for deletion in confirmed_deletions:
        w.writerow([deletion[0], deletion[1], deletion[2], deletion[3], deletion[4]])