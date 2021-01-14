import csv
import sys

file_name = str(sys.argv[1])
output_number = int(sys.argv[2])

with open('split_line_fix.'+output_number, 'w', newline = '') as confirmed:
    w = csv.writer(confirmed, delimiter = '\t')
    if output_number == 1:
        pass
    else output_number != 1:
        with open(file_name, newline = '') as file:
            file_reader = csv.reader(file, delimiter = '\t')
            for row in file_reader:
                current_line = row
                first_read = current_line[3]
                while current_line[3] == current_read:
                    w.writerow(current_line)
                    current_line = next(file_reader)
                break