import csv
import collections
import sys

ltr_file = str(sys.argv[1])

srs_file = str(sys.argv[2])

ltr_and_ltr_output_file = str(sys.argv[3])

ltr_and_other_output_file = str(sys.argv[4])

with open(ltr_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    ltr_srs = [[str(row[0]), int(row[1]), int(row[2]), str(row[3]), int(row[4]), str(row[5])] for row in file_reader]
    
with open(ltr_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    ltr_sr_names = [row[3] for row in file_reader]
    
with open(srs_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    srs_dict = {}
    for row in file_reader:
        if row[3] not in srs_dict:
            srs_dict[row[3]] = [[str(row[0]), int(row[1]), int(row[2]), str(row[3]), int(row[4]), str(row[5])]]
        else:
            srs_dict[row[3]].append([str(row[0]), int(row[1]), int(row[2]), str(row[3]), int(row[4]), str(row[5])])


counted = collections.Counter(ltr_sr_names)
ltr_and_ltr = []
for i in range(len(ltr_sr_names)):
    if counted[ltr_sr_names[i]] == 2 and ltr_sr_names[i] not in ltr_and_ltr:
        ltr_and_ltr.append(ltr_sr_names[i])
ltr_and_ltr_read_locs = []
for i in range(len(ltr_and_ltr)):
    for k in range(len(srs_dict[ltr_and_ltr[i]])):
        ltr_and_ltr_read_locs.append(srs_dict[ltr_and_ltr[i]][k])

with open(ltr_and_ltr_output_file, 'w', newline = '') as output:
    w = csv.writer(output, delimiter = '\t')
    output_list = ltr_and_ltr_read_locs
    for i in range(len(output_list)):
        row = [output_list[i][0], output_list[i][1], output_list[i][2], output_list[i][3], output_list[i][4], output_list[i][5]]
        w.writerow(row)

ltr_and_other = []
for i in range(len(ltr_sr_names)):
    if counted[ltr_sr_names[i]] == 1:
        ltr_and_other.append(ltr_sr_names[i])
ltr_and_other_read_locs = []
for i in range(len(ltr_and_other)):
    for k in range(len(srs_dict[ltr_and_other[i]])):
        ltr_and_other_read_locs.append(srs_dict[ltr_and_other[i]][k])
ltr_and_other_read_locs_tuples = [tuple(l) for l in ltr_and_other_read_locs]
ltr_srs_tuples = [tuple(l) for l in ltr_srs]
ltrs_other_readloc_only = list(set(ltr_srs_and_other_half_tuples) - set(ltr_srs_tuples))

with open(ltr_and_other_output_file, 'w', newline = '') as output:
    w = csv.writer(output, delimiter = '\t')
    output_list = ltrs_other_readloc_only
    for i in range(len(output_list)):
        row = [output_list[i][0], output_list[i][1], output_list[i][2], output_list[i][3], output_list[i][4], output_list[i][5]]
        w.writerow(row)

