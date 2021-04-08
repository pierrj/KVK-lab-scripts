import csv
import sys
import more_itertools as mit
from itertools import filterfalse

def find_ranges(iterable):
    """Yield range of consecutive numbers."""
    for group in mit.consecutive_groups(iterable):
        group = list(group)
        if len(group) == 1:
            yield group[0]
        else:
            yield group[0], group[-1]

deletion_processed_coords = str(sys.argv[1])
gene_loc = str(sys.argv[2])
output_deletion = str(sys.argv[3])
tolerance = int(sys.argv[4])
gene = str(sys.argv[5])
genome = str(sys.argv[6])

with open(deletion_processed_coords, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    deletion_processed_coords_list = []
    for row in file_reader:
        deletion_processed_coords_list.append([row[4], int(row[0]), int(row[1]), row[5], int(row[2]), int(row[3])])

with open(gene_loc, newline = '') as file:
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
for deletion in deletions:
    if deletion[0] <= gene_loc_list[1] and deletion[1] >= gene_loc_list[2]:
        matching_deletion.append(deletion)
if len(matching_deletion) == 0 or len(matching_deletion) > 1:
    raise InputError('matching deletion not found')
single_matching_deletion = matching_deletion[0]

start_corresponding_deletion = 0
end_corresponding_deletion = 0
for match in deletion_processed_coords_list:
    if match[2] == single_matching_deletion[0]-1:
        if match[4] - match[5] > 0:
            start_corresponding_deletion = min([match[4], match[5]])
        elif match[4] - match[5] < 0:
            start_corresponding_deletion = max([match[4], match[5]])
    if match[1] == single_matching_deletion[1]+1:
        if match[4] - match[5] > 0:
            end_corresponding_deletion = max([match[4], match[5]])
        elif match[4] - match[5] < 0:
            end_corresponding_deletion = min([match[4], match[5]])
if abs(end_corresponding_deletion - start_corresponding_deletion) <= tolerance:
    with open(output_deletion , 'w', newline = '') as output_csv:
            w = csv.writer(output_csv, delimiter = '\t') ## lineterminator="\n" might fix issues with windows
            w.writerow([gene_loc_list[0], single_matching_deletion[0], single_matching_deletion[1], gene, genome])