import sys
import csv

test_file = str(sys.argv[1])

with open(test_file, newline = '') as file:
    for row in file_reader:
        print(row)