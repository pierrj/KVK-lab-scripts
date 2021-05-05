## USAGE ##
# this python script uses hierarchical clustering to merge highly similar eccDNA forming regions together
# input format is currently hard coded and should be tab delimited file that looks like so:
## first row is header
## columns are Triplet, Amino Acid, Blank, Blank, Blank, Triplet, Fraction(yeast), Frequency(yeast), Number(yeast), Fraction(ecoli), Frequency(ecoli), Number(ecoli), Fraction(rice), Frequency(rice), Number(rice), Fraction(benthi), Frequency(benthi), Number(benthi)
# inputs (in order)
## input table
## output name
## minimum fraction codon usage of other organism (not yeast) i.e. 0.1
## maximum fraction to remove of codon usage of yeast i.e. 0.3
#
# Example usage
## python3 ./rebalance_codon_usage.py codon_usage.tsv codon_usage_rebalanced.tsv 0.1 0.3


import csv
import numpy as np
import sys

codon_table = sys.argv[1]
output = sys.argv[2]
min_fraction_any = int(sys.argv[3])
max_remove_fraction = int(sys.argv[4])

## read codon table into dictionary of list of lists, one list of lists for each amino acid
codon_table_dict = {}
with open(codon_table, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    next(file_reader) ## skip header
    for row in file_reader:
        if row[1] not in codon_table_dict:
            codon_table_dict[row[1]] = []
        codon_table_dict[row[1]].append([row[0], float(row[6]), float(row[7]), int(row[8]), 
                                         float(row[10]), float(row[11]), int(row[12]),
                                        float(row[14]), float(row[15]), int(row[16]),
                                        float(row[18]), float(row[19]), int(row[20])])

## convert to numpy 2D array
codon_table_arrays = {}
for key in codon_table_dict.keys():
    codon_table_arrays[key] = np.array(codon_table_dict[key], dtype=object)

codon_table_optimized = {}
for key in codon_table_arrays:
    array = codon_table_arrays[key]
    ## fractions in all three other organisms must be greater than min_fraction_any
    ## or fraction in first organism is greater than max_remove_fraction
    subset_array = array[np.logical_or(np.logical_and.reduce((array[:,4] >= min_fraction_any, array[:,7] >= min_fraction_any, array[:,10] >= min_fraction_any)),
                                       array[:,1] >= max_remove_fraction)]
    ## make sure there is at least one left?
    if np.shape(subset_array)[0] == 0:
        print(key + ' is empty after filtering, keeping all codons')
        codon_table_optimized[key] = array
    else:
        codon_table_optimized[key] = subset_array

## recalculate fractions for each aa based off new totals and add to list
codon_table_final = []
for key in codon_table_optimized:
    array = codon_table_optimized[key]
    sum_one = np.sum(array[:,3])
    sum_two = np.sum(array[:,6])
    sum_three = np.sum(array[:,9])
    sum_four = np.sum(array[:,12])
    for line in array:
        codon_table_final.append([line[0], key, '', '', '', line[0], round(line[3]/sum_one,2), line[2] , line[3],
                                 round(line[6]/sum_two,2), line[5] , line[6],
                                 round(line[9]/sum_three,2), line[8] , line[9],
                                 round(line[12]/sum_four,2), line[11], line[12]])

## write list
with open(output, 'w', newline = '') as output_file:
    w = csv.writer(output_file, delimiter = '\t')
    for row in codon_table_final:
        w.writerow(row)