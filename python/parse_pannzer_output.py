#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
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