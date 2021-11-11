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

## USAGE ##
# this script takes the locations of split reads that intersect with the LTR region of a retrotransposon on one side and finds the second side of the split read
# "ltr_file" bed file of ltr locations
# "srs_file" bed file of sr locations
# "ltr_and_ltr_output_file" where reads that intersect LTRs on both sides are printed
# "ltr_and_other_output_file" where reads that intersect LTR and something else on the otehr side are printed

ltr_file = str(sys.argv[1])

srs_file = str(sys.argv[2])

ltr_and_ltr_output_file = str(sys.argv[3])

ltr_and_other_output_file = str(sys.argv[4])

# read in files as lists or dictionaries
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


# simple count to see if the reads overlap ltrs twice
counted = collections.Counter(ltr_sr_names)
ltr_and_ltr = []
for i in range(len(ltr_sr_names)):
    if counted[ltr_sr_names[i]] == 2 and ltr_sr_names[i] not in ltr_and_ltr:
        ltr_and_ltr.append(ltr_sr_names[i])
ltr_and_ltr_read_locs = []
for i in range(len(ltr_and_ltr)):
    for k in range(len(srs_dict[ltr_and_ltr[i]])):
        ltr_and_ltr_read_locs.append(srs_dict[ltr_and_ltr[i]][k])

# write
with open(ltr_and_ltr_output_file, 'w', newline = '') as output:
    w = csv.writer(output, delimiter = '\t')
    output_list = ltr_and_ltr_read_locs
    for i in range(len(output_list)):
        row = [output_list[i][0], output_list[i][1], output_list[i][2], output_list[i][3], output_list[i][4], output_list[i][5]]
        w.writerow(row)

# if they only appear once then that means that the second side of the split read overlaps something else
ltr_and_other = []
for i in range(len(ltr_sr_names)):
    if counted[ltr_sr_names[i]] == 1:
        ltr_and_other.append(ltr_sr_names[i])
ltr_and_other_read_locs = []
for i in range(len(ltr_and_other)):
    for k in range(len(srs_dict[ltr_and_other[i]])):
        ltr_and_other_read_locs.append(srs_dict[ltr_and_other[i]][k])
# tuple stuff because of the dictionary...
ltr_and_other_read_locs_tuples = [tuple(l) for l in ltr_and_other_read_locs]
ltr_srs_tuples = [tuple(l) for l in ltr_srs]
ltrs_other_readloc_only = list(set(ltr_and_other_read_locs_tuples) - set(ltr_srs_tuples))

# write
with open(ltr_and_other_output_file, 'w', newline = '') as output:
    w = csv.writer(output, delimiter = '\t')
    output_list = ltrs_other_readloc_only
    for i in range(len(output_list)):
        row = [output_list[i][0], output_list[i][1], output_list[i][2], output_list[i][3], output_list[i][4], output_list[i][5]]
        w.writerow(row)

