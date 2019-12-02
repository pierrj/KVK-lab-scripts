import csv
import statistics

with open('parallel.confirmed') as confirmed:
    confirmed_reader = csv.reader(confirmed, delimiter = '\t')
    confirmed_list = [[int(row[0]), int(row[1]), int(row[2])] for row in confirmed_reader]

with open('genomecoverage.mergedandpe.G3_1A_bwamem.bed') as coverage:
    coverage_reader = csv.reader(coverage, delimiter = '\t')
    coverage_indexed = [[] for i in range(56)]
    for row in coverage_reader:
        coverage_indexed[(int(row[0][10:12])-1)].append([int(row[1]) -1, int(row[2])])

for i in range(len(confirmed_list)):
    ecc = confirmed_list[i]
    region_len = ecc[2] - ecc[1]
    region = coverage_indexed[ecc[0]][ecc[1]:(ecc[2]+1)]
    region_cov = [region[i][1] for i in range(len(region))]
    beforestart = ecc[1] - 1 - region_len
    if beforestart <= 0:
        ecc.append('lconf')
        ecc.append('too close to start')
        continue
    afterstart = ecc[2] + 1
    if afterstart + region_len > coverage_indexed[ecc[0]][-1][0]:
        ecc.append('lconf')
        ecc.append('too close to end')
        continue
    region_before = coverage_indexed[ecc[0]][beforestart:beforestart + 1 + region_len]
    region_before_cov = [region_before[i][1] for i in range(len(region_before))]
    region_after = coverage_indexed[ecc[0]][afterstart:afterstart + 1 + region_len]
    region_after_cov = [region_after[i][1] for i in range(len(region_after))]
    mean_region = statistics.mean(region_cov)
    mean_before = statistics.mean(region_before_cov)
    mean_after = statistics.mean(region_after_cov)
    if mean_region >= (2*mean_before) and mean_region >= (2*mean_after):
        ecc.append('hconf')
        ecc.append('none')
    else:
        ecc.append('lconf')
        ecc.append('coverage too low')

with open('ecc_confidence', 'w', newline = '') as conf:
    w = csv.writer(conf, delimiter = '\t')
    w.writerows(confirmed_list)