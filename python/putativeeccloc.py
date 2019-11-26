import csv
import sys
with open(str(sys.argv[1]), newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    with open(str(sys.argv[2]), 'w', newline = '') as filtered:
        w = csv.writer(filtered, delimiter = '\t')
        for line1 in file_reader:
            line2 = next(file_reader)
            line1[2] = line2[1]
            w.writerow(line1)
