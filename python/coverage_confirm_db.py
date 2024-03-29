import sys
import csv
import sqlite3
import statistics

output_name = str(sys.argv[1])

output_number = str(sys.argv[2])

with open('merged.confirmed'+output_number) as merged:
    merged_reader = csv.reader(merged, delimiter = '\t')
    flat_merged_list = [[int(row[0]), int(row[1]), int(row[2]), int(row[3]), str(row[4])] for row in merged_reader]


def confidence_check(ecc):
    # get coverage of confirmed ecc region
    region_len = ecc[2] - ecc[1]
    beforestart = ecc[1] - region_len
    afterstart =  ecc[2] + 2 + region_len
    if beforestart > 0:
        region = [ecc[0], beforestart, afterstart]
    else:
        region = [ecc[0], ecc[1], afterstart]
    conn = sqlite3.connect(r"scaffold"+str(region[0]+1)+"_sql.db")
    c = conn.cursor()
    c.execute("SELECT * FROM server WHERE base BETWEEN "+str(region[1])+" AND " +str(region[2]))
    region_all = c.fetchall()
    conn.close()
    region_cov = [region_all[k][1] for k in range(len(region_all))]
    if beforestart > 0:
        ecc_region_cov = region_cov[region_len+1:((2 * region_len+1)+1)]
        before_region_cov = region_cov[0:region_len+1]
        after_region_cov = region_cov[(2 * region_len)+2:]
    else:
        ecc_region_cov = region_cov[0:region_len+1]
        after_region_cov = region_cov[region_len+1:((2 * region_len+1)+1)]
    mean_region = round(statistics.mean(ecc_region_cov), 2)
    if beforestart > 0:
        mean_before = round(statistics.mean(before_region_cov), 2)
    else:
        mean_before = 'N/A'
    if len(after_region_cov) == region_len + 1:
        mean_after = round(statistics.mean(after_region_cov), 2)
    else:
        mean_after = 'N/A'
    # write coverage string containing the means of the regions before, within, and after the confirmed ecc
    coverage_string = str(mean_before)+';'+str(mean_region)+';'+str(mean_after)
    # if less than 95% of the ecc region is covered than the ecc is low confidence
    if ecc_region_cov.count(0) / len(ecc_region_cov) > 0.05:
        ecc.append('lowq')
        ecc.append('incomplete_coverage')
        ecc.append(coverage_string)
        return ecc
    if ecc[3] == 1:
        ecc.append('lowq')
        ecc.append('only_one_splitread')
        ecc.append(coverage_string)
        return ecc
    # if the ecc before or after regions fall beyond the borders of the chromosome than the ecc is medium confidence
    if beforestart <= 0 and ecc[3] < 3:
        ecc.append('conf')
        ecc.append('too_close_to_start')
        ecc.append(coverage_string)
        return ecc
    if len(after_region_cov) != region_len + 1 and ecc[3] <3:
        ecc.append('conf')
        ecc.append('too_close_to_end')
        ecc.append(coverage_string)
        return ecc
    if beforestart <= 0 and ecc[3] >= 3:
        ecc.append('hconf')
        ecc.append('splitreads')
        ecc.append(coverage_string)
        return ecc
    if len(after_region_cov) != region_len + 1 and ecc[3] >= 3:
        ecc.append('hconf')
        ecc.append('splitreads')
        ecc.append(coverage_string)
        return ecc
    # if the mean coverage of the eccDNA is twice the size of the before and after regions OR the eccDNA is supported by three or more splitreads then eccDNA is high confidence, otherwise is it medium confidence
    if mean_region >= (2*mean_before) and mean_region >= (2*mean_after) and ecc[3] < 3:
        ecc.append('hconf')
        ecc.append('coverage')
        ecc.append(coverage_string)
        return ecc
    if (mean_region < (2*mean_before) or mean_region < (2*mean_after)) and ecc[3] >= 3:
        ecc.append('hconf')
        ecc.append('splitreads')
        ecc.append(coverage_string)
        return ecc
    if mean_region >= (2*mean_before) and mean_region >= (2*mean_after) and ecc[3] >= 3:
        ecc.append('hconf')
        ecc.append('splitreads_and_coverage')
        ecc.append(coverage_string)
    else:
        ecc.append('conf')
        ecc.append('coverage_and_splitreads_too_low')
        ecc.append(coverage_string)
    return ecc

confidence_flat_merged_list = list(map(confidence_check, flat_merged_list))

with open('ecccaller_output.' + output_name + '.details.tsv'+output_number, 'w', newline = '') as bed:
    w = csv.writer(bed, delimiter = '\t')
    for i in range(len(confidence_flat_merged_list)):
        row = [confidence_flat_merged_list[i][0]+1, confidence_flat_merged_list[i][1], confidence_flat_merged_list[i][2], confidence_flat_merged_list[i][3], confidence_flat_merged_list[i][4],confidence_flat_merged_list[i][5], confidence_flat_merged_list[i][6], confidence_flat_merged_list[i][7]]
        w.writerow(row)

# write bed file output, with actual chromosome names, start and end coordinates and ecc number
with open('ecccaller_output.' + output_name + '.bed'+output_number, 'w', newline = '') as bed:
    w = csv.writer(bed, delimiter = '\t')
    for i in range(len(confidence_flat_merged_list)):
        if confidence_flat_merged_list[i][5] == 'lowq':
            color = '255,0,0'
        if confidence_flat_merged_list[i][5] == 'conf':
            color = '255,255,0'
        if confidence_flat_merged_list[i][5] == 'hconf':
            color = '0,255,0'
        row = [confidence_flat_merged_list[i][0]+1, confidence_flat_merged_list[i][1], confidence_flat_merged_list[i][2], i, 0, '+', confidence_flat_merged_list[i][1], confidence_flat_merged_list[i][2], color ]
        w.writerow(row)