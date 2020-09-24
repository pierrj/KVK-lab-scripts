import re as regex
import csv
import numpy as np
import collections
import sys


def regex_partition(content, separator):
    separator_match = regex.search(separator, content)
    if not separator_match:
        return content, content, content

    matched_separator = separator_match.group(0)
    parts = regex.split(matched_separator, content, 1)

    return parts[0], matched_separator, parts[1]

def trim_domain_names(input_list):
    regexp = '_[0-9]'
    trimmed_domains = []
    for i in range(len(input_list)):
        if len(regex_partition(input_list[i], regexp)[0]) > 2:
            trimmed_domains.append(regex_partition(input_list[i], regexp)[0])
        else:
            trimmed_domains.append(input_list[i])
    regexp = '_C|_N'
    trimmed_domains_noCN = []
    for i in range(len(trimmed_domains)):
        trimmed_domains_noCN.append(regex_partition(trimmed_domains[i], regexp)[0])
        if len(regex_partition(trimmed_domains[i], regexp)[0]) > 2:
            trimmed_domains_noCN.append(regex_partition(trimmed_domains[i], regexp)[0])
        else:
            trimmed_domains_noCN.append(trimmed_domains[i])
    regexp = 'Peptidase'
    trimmed_domains_noCN_peptidasegrouped = []
    for i in range(len(trimmed_domains_noCN)):
        trimmed_domains_noCN_peptidasegrouped.append(regex_partition(trimmed_domains_noCN[i], regexp)[1])
    return trimmed_domains_noCN_peptidasegrouped

def count(input_list):
    counted = collections.Counter(input_list)
    counted_unique = [[key, counted[key]] for key in counted]
    return counted_unique

gene_subset = str(sys.argv[1])
all_genes = str(sys.argv[2])

with open(gene_subset, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    G3_list = []
    for row in file_reader:
        split = row[1].split("~")
        for domain in split:
            G3_list.append(domain)

with open(all_genes, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    WG_list = []
    for row in file_reader:
        split = row[1].split("~")
        for domain in split:
            WG_list.append(domain)

trimmed_G3 = trim_domain_names(G3_list)
counted_G3 = collections.Counter(trimmed_G3)
counted_trimmed_G3 = [[key, counted_G3[key]] for key in counted_G3]
trimmed_WG = trim_domain_names(WG_list)
counted_WG = collections.Counter(trimmed_WG)
WG_dict = {key:counted_WG[key] for key in counted_WG}

normalized_G3 = []
for i in range(len(counted_trimmed_G3)):
    domain = counted_trimmed_G3[i][0]
    g3_count = counted_trimmed_G3[i][1]
    wg_count = WG_dict[domain]
    if wg_count >= 10:
        fraction = round(g3_count/wg_count*100)
        normalized_G3.append([domain, fraction])
sorted_normalized_G3 = sorted(normalized_G3, reverse=True, key=lambda x: x[1])

no_normalized_G3 = []
for i in range(len(counted_trimmed_G3)):
    domain = counted_trimmed_G3[i][0]
    g3_count = counted_trimmed_G3[i][1]
    wg_count = WG_dict[domain]
    if wg_count >= 10:
        no_normalized_G3.append([domain, g3_count])
sorted_no_normalized_G3 = sorted(no_normalized_G3, reverse=True, key=lambda x: x[1])

total_output = str(sys.argv[3])

normalized_output = str(sys.argv[4])


with open(total_output, 'w', newline = '') as file:
    w = csv.writer(file, delimiter = '\t')
    for i in range(len(sorted_no_normalized_G3)):
        row = [sorted_no_normalized_G3[i][0], sorted_no_normalized_G3[i][1]]
        w.writerow(row)
        
with open(normalized_output, 'w', newline = '') as file:
    w = csv.writer(file, delimiter = '\t')
    for i in range(len(sorted_normalized_G3)):
        row = [sorted_normalized_G3[i][0], sorted_normalized_G3[i][1]]
        w.writerow(row)