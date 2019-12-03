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
    region_cov = [region[k][1] for k in range(len(region))]
    mean_region = round(statistics.mean(region_cov), 2)
    beforestart = ecc[1] - 1 - region_len
    afterstart = ecc[2] + 1
    region_before = coverage_indexed[ecc[0]][beforestart:beforestart + 1 + region_len]
    region_before_cov = [region_before[j][1] for j in range(len(region_before))]
    region_after = coverage_indexed[ecc[0]][afterstart:afterstart + 1 + region_len]
    region_after_cov = [region_after[g][1] for g in range(len(region_after))]
    if beforestart > 0:
        mean_before = round(statistics.mean(region_before_cov), 2)
    else:
        mean_before = 'N/A'
    if afterstart + region_len <= coverage_indexed[ecc[0]][-1][0]:
        mean_after = round(statistics.mean(region_after_cov), 2)
    else:
        mean_after = 'N/A'
    coverage_string = str(mean_before)+';'+str(mean_region)+';'+str(mean_after)
    if region_cov.count(0) / len(region_cov) > 0.05:
        ecc.append('lowq')
        ecc.append('incomplete_coverage')
        ecc.append(coverage_string)
        continue
    if beforestart <= 0:
        ecc.append('conf')
        ecc.append('too_close_to_start')
        ecc.append(coverage_string)
        continue
    if afterstart + region_len > coverage_indexed[ecc[0]][-1][0]:
        ecc.append('conf')
        ecc.append('too_close_to_end')
        ecc.append(coverage_string)
        continue
    if mean_region >= (2*mean_before) and mean_region >= (2*mean_after):
        ecc.append('hconf')
        ecc.append('none')
        ecc.append(coverage_string)
    else:
        ecc.append('conf')
        ecc.append('coverage_too_low')
        ecc.append(coverage_string)

with open('ecc_confidence', 'w', newline = '') as conf:
    w = csv.writer(conf, delimiter = '\t')
    w.writerows(confirmed_list)