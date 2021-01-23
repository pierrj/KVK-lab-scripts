import csv
import re
import random
from scipy.stats import gaussian_kde
import numpy as np
import statsmodels.api as sm
import scipy.interpolate
import sys

known_srs_length = str(sys.argv[1])
split_reads = str(sys.argv[2])
column_cutoff = int(sys.argv[3])
output_number = str(sys.argv[4])

with open(known_srs_length, newline = '') as file: ## this is for each individual technical replicate
    file_reader = csv.reader(file, delimiter = '\t')
    dsn_list = [int(row[0]) for row in file_reader]
density = gaussian_kde(dsn_list)
kde_max = 0
while True:
    if density(kde_max) == 0:
        break
    else:
        kde_max += 100
xs = np.linspace(0,kde_max,1000)
densities = density(xs)
y_interp = scipy.interpolate.interp1d(xs,densities)

## define process split reads
start_pattern = "^([0-9]+)M.*[HS]$"
end_pattern = ".*[HS]([0-9]+)M$"

def process_split_read(read):
    matches = re.findall(r'(\d+)([A-Z]{1})', read[6])
    matches_sums = {'M': 0, 'other': 0}
    for i in range(len(matches)):
        if matches[i][1] == 'M':
            matches_sums['M'] += int(matches[i][0])
        else:
            matches_sums['other'] += int(matches[i][0])
    if str(read[5]) == '-':
        sense = 'reverse'
    elif str(read[5]) == '+':
        sense = 'forward'
    string = read[6]
    if sense == 'forward':
        if re.match(start_pattern, string):
            loc = 'start'
        elif re.match(end_pattern, string):
            loc = 'end'
        else:
            return False
    elif sense == 'reverse':
        if re.match(start_pattern, string):
            loc = 'end'
        elif re.match(end_pattern, string):
            loc = 'start'
        else:
            return False ## if regexp doesn't work then just drop the split
    if matches_sums['M'] >= matches_sums['other']:
        split_read_side_one.append([read[0], int(read[1]), int(read[2]), sense, read[6], loc, matches_sums['M'], matches_sums['other']])
    elif matches_sums['M'] < matches_sums['other']:
        split_read_side_two.append([read[0], int(read[1]), int(read[2]), sense, read[6], loc, matches_sums['M'], matches_sums['other']])
    return True

def choose_split_reads(split_read_side_one, split_read_side_two_pre):
    if len(split_read_side_two_pre) == 0:
        return False
    split_read_side_two = np.array(split_read_side_two_pre, dtype=object)
    combos = []
    for i in range(len(split_read_side_one)):
        split_read_first_side = split_read_side_one[i]
        ## reduce to
        # same chromosome
        # same orientation
        # distance smaller than column cutoff
        # start and end pairs
        # properly oriented (eccDNA vs intron)
        # total of read lengths is around 150
        if ((split_read_first_side[3] == 'forward' and split_read_first_side[5] == 'start') or 
            (split_read_first_side[3] == 'reverse' and split_read_first_side[5] == 'end')):
                reduced = split_read_side_two[np.logical_and.reduce((split_read_side_two[:, 0] == split_read_first_side[0],
                                       split_read_side_two[:, 3] == split_read_first_side[3],
                                       abs(split_read_side_two[:,1] - split_read_first_side[1]) < column_cutoff,
                                        split_read_side_two[:, 5] != split_read_first_side[5],
                                        split_read_side_two[:,1] < split_read_first_side[1],
                                        np.isclose((split_read_side_two[:,6] + split_read_first_side[6]).astype(int), (split_read_side_two[:,7] + split_read_first_side[7]).astype(int), atol=10))]
        if ((split_read_first_side[3] == 'forward' and split_read_first_side[5] == 'end') or
            (split_read_first_side[3] == 'reverse' and split_read_first_side[5] == 'start')):
                reduced = split_read_side_two[np.logical_and.reduce((split_read_side_two[:, 0] == split_read_first_side[0],
                                       split_read_side_two[:, 3] == split_read_first_side[3],
                                       abs(split_read_side_two[:,1] - split_read_first_side[1]) < column_cutoff,
                                        split_read_side_two[:, 5] != split_read_first_side[5],
                                        split_read_side_two[:,1] > split_read_first_side[1],
                                        np.isclose((split_read_side_two[:,6] + split_read_first_side[6]).astype(int), (split_read_side_two[:,7] + split_read_first_side[7]).astype(int), atol=10)))]
        if len(reduced) != 0:
            for ii in range(len(reduced)):
                split_read_second_side = reduced[ii]
                start = min([split_read_first_side[1], split_read_first_side[2], split_read_second_side[1], split_read_second_side[2]])
                end = max([split_read_first_side[1], split_read_first_side[2], split_read_second_side[1], split_read_second_side[2]])
                combos.append([split_read_first_side[0], start, end])
    if combos == []: # any of the options left?
        return False
    densities = []
    for ii in range(len(combos)):
        try:
            densities.append(y_interp(abs(combos[ii][2]-combos[ii][1])))
        except ValueError:
            densities.append(0)
    densities_sum = sum(densities)
    if densities_sum == 0: ## any of the options possible?
        return False
    densities_array = np.array(densities)
    densities_ratio = densities_array/densities_sum
    choice = random.choices(combos, densities_ratio, k=1)[0]
    return choice

with open(split_reads, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    with open('mapq0_choices.'+output_number, 'w', newline = '') as confirmed:
        w = csv.writer(confirmed, delimiter = '\t')
        for row in file_reader:
            current_line = row
            current_read = current_line[3]
            split_read_side_one = []
            split_read_side_two = []
            while current_line[3] == current_read:
                process_split_read(current_line)
                try:
                    current_line = next(file_reader)
                except StopIteration:
                    break
            choice = choose_split_reads(split_read_side_one, split_read_side_two)
            if choice:
                w.writerow([choice[0], choice[1], choice[2]])