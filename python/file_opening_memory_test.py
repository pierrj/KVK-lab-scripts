import csv
import ipyparallel as ipp
import sys

split_read_file = str(sys.argv[1])

outwardfacing_read_file = str(sys.argv[2])

output_name = str(sys.argv[3])

scaffold_number = int(sys.argv[4])

print('successfully load modules')

# open putative ecc list and index to speed up confirming eccs
with open(split_read_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    eccloc_list = []
    for row in file_reader:
        ecc_loc = [int(row[0]) - 1, int(row[1]), int(row[2])]
        eccloc_list.append(ecc_loc)

print('successfully opened split read file')

# open opposite facing discordant read file
with open(outwardfacing_read_file, newline = '') as discordant:
    discordant_reader = csv.reader(discordant, delimiter = '\t')
    # index discordant read file so that confirmeccs() only looks at discordant reads found on the same chromosome
    discordant_indexed = [[] for i in range(scaffold_number)]
    for row in discordant_reader:
        discordant_indexed[(int(row[0])-1)].append([int(row[1]), int(row[2])])

print('successfully opened files')