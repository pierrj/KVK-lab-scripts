import csv
import sys
import math
import collections

column_2 = str(sys.argv[1])
bin_width = int(sys.argv[2])
max_bin = int(sys.argv[3])
output = str(sys.argv[4])

with open(column_2, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    column_2_list = []
    for row in file_reader:
        column_2_list.append([int(row[0])])

column_2_round = []
for i in range(len(column_2_list)):
    rounded = int(math.floor(column_2_list[i][0] / 100.0)) * 100
    column_2_round.append(rounded)

counted = collections.Counter(column_2_round)
count_list = []
for x in set(column_2_round):
    count_list.append([x, counted[x]])
count_list_sorted = sorted(count_list, key=lambda x: x[0])

while count_list_sorted[-1][0] != max_bin:
    count_list_sorted.append([count_list_sorted[-1][0]+100, 0])

total_eccs = 0
for i in range(len(count_list_sorted)):
    total_eccs += count_list_sorted[i][1]
count_list_sorted_density = []
for i in range(len(count_list_sorted)):
    density = count_list_sorted[i][1]/total_eccs
    count_list_sorted_density.append([count_list_sorted[i][0], density])

with open(output, 'w', newline = '') as output_csv:
    w = csv.writer(output_csv, delimiter = '\t')
    for i in range(len(count_list_sorted_density)):
        row = [count_list_sorted_density[i][0], count_list_sorted_density[i][1]]
        w.writerow(row)