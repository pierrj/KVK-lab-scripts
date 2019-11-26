import csv
import sys

with open(str(sys.argv[1]), newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    with open(str(sys.argv[2]), 'w', newline = '') as filtered:
        w = csv.writer(filtered, delimiter = '\t')
        for line1 in file_reader:
            line2 = next(file_reader)
            sum = int(line1[3]) - int(line2[3])
            if abs(sum) <= 1000000:
                w.writerow(line1)
                w.writerow(line2)
