import csv
import collections
import sys

input_file = str(sys.argv[1])
score_cutoff = float(sys.argv[2])
output_name = str(sys.argv[3])

with open(input_file) as file:
    file_reader = csv.reader(file, delimiter = '\t')
    pannzer_output = collections.defaultdict(list)
    for row in file_reader:
        if row and row[0] != 'qpid' and float(row[5]) > score_cutoff:
            pannzer_output[row[0]].append('GO:' + row[2])

with open(output_name, 'w', newline='') as csv_file:  
    writer = csv.writer(csv_file, delimiter = '\t')
    for key, value in pannzer_output.items():
        if len(value) > 1:
            to_write = [key]
            to_join = []
            for i in range(len(value)):
                to_join.append(value[i]+",")
            to_write.append(" ".join(to_join)[:-1])
        else:
            to_write = [key]
            to_write.append(value[0])
        writer.writerow(to_write)