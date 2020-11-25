import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
import csv
import sys

outputname = str(sys.argv[1])
scaffold_number = int(sys.argv[2])
max_d = int(sys.argv[3])

def get_distance_from_point(row):
    start_distance = abs(row['start'] - start_mean)
    end_distance = abs(row['end'] - end_mean)
    return start_distance + end_distance

parallel_confirmed = pd.read_csv("parallel.confirmed", sep = '\t', names = ['chrom', 'start', 'end'])

parallel_confirmed = parallel_confirmed.groupby(parallel_confirmed.columns.tolist()).size().reset_index().\
    rename(columns={0:'splitreads'})

parallel_confirmed_dict = {}
for key in range(scaffold_number):
    parallel_confirmed_dict[key] = parallel_confirmed.loc[parallel_confirmed['chrom'] == key]

representative_eccs= []
representative_variants = {}
for chrom in range(scaffold_number):
    scaffold_subset = parallel_confirmed_dict[chrom]
    scaffold_subset_startend =  scaffold_subset[["start", "end"]]
    if scaffold_subset_startend.empty:
        continue
    else:
        linkage_for_scaffold = linkage(scaffold_subset_startend) ## need to verify clustering method to use here
    clusters = fcluster(linkage_for_scaffold, max_d, criterion='distance')
    for i in range(1, np.amax(clusters)+1):
        boolean_array = clusters==i
        cluster_subset = scaffold_subset[boolean_array] ## this should be the current scaffold
        if len(cluster_subset) == 1:
            clustered = 'no'
            rep_start = cluster_subset.iloc[0]['start']
            rep_end = cluster_subset.iloc[0]['end']
            split_read_count = cluster_subset.iloc[0]['splitreads']
            point_of_interest = [chrom, rep_start, rep_end, split_read_count ,clustered]
            representative_eccs.append(point_of_interest)
        else:
            clustered = 'yes'
            start_mean = cluster_subset.mean(axis=0)['start']
            end_mean = cluster_subset.mean(axis=0)['end']
            cluster_subset['distance_from_center'] = cluster_subset.apply(get_distance_from_point, axis=1)
            rep_start = int(cluster_subset[cluster_subset.distance_from_center == cluster_subset.distance_from_center.min()].iloc[0]['start'])
            rep_end = int(cluster_subset[cluster_subset.distance_from_center == cluster_subset.distance_from_center.min()].iloc[0]['end'])
            split_read_count = cluster_subset['splitreads'].sum()
            point_of_interest = [chrom, rep_start, rep_end, split_read_count ,clustered]
            representative_eccs.append(point_of_interest)
            cluster_subset = cluster_subset[["chrom","start", "end", "splitreads"]]
            cluster_subset_listoflists = cluster_subset.values.tolist()
            cluster_subset_tupleoftuples = tuple(tuple(sub) for sub in cluster_subset_listoflists)
            representative_variants[tuple(point_of_interest)] = cluster_subset_tupleoftuples

with open('merged.confirmed'+outputname, 'w', newline = '') as merged:
    w = csv.writer(merged, delimiter = '\t')
    w.writerows(representative_eccs)

with open('ecccaller_splitreads.' + outputname + '.tsv', 'w', newline="") as variants_dict:
    w = csv.writer(variants_dict, delimiter = '\t')
    for key, value in representative_variants.items():
        w.writerow([key, value])