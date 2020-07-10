import sys
import csv

test_file = str(sys.argv[1])

with open(test_file, newline = '') as file:
    list1 = []
    for row in file_reader:
        list1.append(row)

output = 'processed' + test_file

with open(output, 'w', newline = '') as out:
    w = csv.writer(out, delimiter = '\t')
    w.writerows(out)