import sys
import csv
import re
import math

input_file = str(sys.argv[1])
output_file = str(sys.argv[2])

with open(input_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    with open(output_file, 'w', newline = '') as confirmed:
        w = csv.writer(confirmed, delimiter = '\t')
        for row1 in file_reader:
            row2 = next(file_reader)
            matches_row1 = re.findall(r'(\d+)([A-Z]{1})', row1[5])
            matches_sums_row1 = {'M': 0, 'other': 0}
            for i in range(len(matches_row1)):
                if matches_row1[i][1] == 'M':
                    matches_sums_row1['M'] += int(matches_row1[i][0])
                else:
                    matches_sums_row1['other'] += int(matches_row1[i][0])
            matches_row2 = re.findall(r'(\d+)([A-Z]{1})', row2[5])
            matches_sums_row2 = {'M': 0, 'other': 0}
            for i in range(len(matches_sums_row2)):
                if matches_row2[i][1] == 'M':
                    matches_sums_row2['M'] += int(matches_row2[i][0])
                else:
                    matches_sums_row2['other'] += int(matches_row2[i][0])
            read_length_matches = matches_sums_row1['M'] + matches_sums_row2['M']
            read_length_nonmatches = matches_sums_row1['other'] + matches_sums_row2['other']
            if math.isclose(read_length_matches,read_length_nonmatches, abs_tol=10):
                w.writerow(row1)
                w.writerow(row2)