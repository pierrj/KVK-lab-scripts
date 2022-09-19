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
import numpy as np
import sys

subdir = sys.argv[1]

all_samples = ['G3_1A', 'G3_1B', 'G3_1C', 'G3_2A', 'G3_2B', 'G3_2C', 'G3_3A', 'G3_3B']
biorep_1 = ['G3_1A', 'G3_1B', 'G3_1C']
biorep_2 = ['G3_2A', 'G3_2B', 'G3_2C']
biorep_3 = ['G3_3A', 'G3_3B']
all_bioreps = [biorep_1, biorep_2, biorep_3]

## make dictionaries of eccdna forming regions
ecc_dict = {}

for sample in all_samples:
    # print(sample)
    eccs_bed = []
    with open(subdir+'/'+sample+'.ecc_caller_out.splitreads.bed', newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        for row in file_reader:
            eccs_bed.append([row[0], row[1], row[2]])
    eccs_bed_indexed = [[] for i in range(56)]
    for ecc in eccs_bed:
        scaffold_num = int(ecc[0][10:12])-1
        eccs_bed_indexed[scaffold_num].append([int(ecc[1]), int(ecc[2])])
    eccs_bed_arrays = []
    for i in range(len(eccs_bed_indexed)):
        arrays = np.array(eccs_bed_indexed[i])
        if arrays.size != 0: # some scaffolds have no eccs
            eccs_bed_arrays.append(np.unique(arrays,axis=0)) # unique gets rid of split reads here
        else:
            eccs_bed_arrays.append(arrays)
    ecc_dict[sample] = eccs_bed_arrays

list_of_lists_output = []

for tolerance in [0, 10, 50, 100, 250, 500, 1000]: ## check different tolerances
    # print('for tolerance')
    # print(tolerance)
    final_dict={}
    # print('getting data for bio reps...')
    for biorep in all_bioreps:
        biorep_string = biorep[0][0:4]
        # print(biorep_string)
        final_dict[biorep_string] = []
        for i in range(56):
            samples_list = []
            for sample in biorep:
                if ecc_dict[sample][i].size != 0:
                    samples_list.append(ecc_dict[sample][i])
            if len(samples_list) != 0:
                samples_concatenated = np.concatenate((samples_list), axis=0) ## concatenate all eccs in all tech reps
            if samples_concatenated.size != 0:
                unique_samples_concatenated = np.unique(samples_concatenated, axis=0)
            else:
                unique_samples_concatenated = samples_concatenated
            if unique_samples_concatenated.size != 0:
                for sample in biorep:
                    sample_eccs = ecc_dict[sample][i]
                    membership = []
                    for ecc in unique_samples_concatenated: # loop through all eccs
                        # check to see if any ecc overlaps in that tech rep
                        if sample_eccs.size != 0:
                            if np.size(
                            sample_eccs[np.logical_and(
                                np.isclose(sample_eccs[:,0], ecc[0], atol=tolerance, rtol=0),
                                np.isclose(sample_eccs[:,1], ecc[1], atol=tolerance, rtol=0)
                            )] ## add isclose or just == here, depending on the level of tolerance
                            ) != 0:
                                membership.append(1)
                            else:
                                membership.append(0)
                        else:
                            membership.append(0)
                    membership_column = np.reshape(np.array(membership), (-1, 1))
                    unique_samples_concatenated = np.append(unique_samples_concatenated, membership_column, axis=1)
            final_dict[biorep_string].append(unique_samples_concatenated)
    ## counts for overlap between tech reps
    # print('overlaps per biorep are...')
    # print('(order is A, B, C, AB, BC, AC, ABC)')
    for biorep in final_dict.keys():
        A_count = 0
        B_count = 0
        C_count = 0
        AB_count = 0
        BC_count = 0
        AC_count = 0
        ABC_count = 0
        for i in range(56):
            eccs = final_dict[biorep][i]
            if eccs.size != 0:
                if biorep != 'G3_3':
                    A_count += len(eccs[(eccs[:,2] == 1) & 
                                        (eccs[:,3] == 0) & 
                                        (eccs[:,4] == 0)
                    ])
                    B_count += len(eccs[(eccs[:,2] == 0) & 
                                        (eccs[:,3] == 1) & 
                                        (eccs[:,4] == 0)
                    ])
                    C_count += len(eccs[(eccs[:,2] == 0) & 
                                        (eccs[:,3] == 0) & 
                                        (eccs[:,4] == 1)
                    ])
                    AB_count += len(eccs[(eccs[:,2] == 1) & 
                                        (eccs[:,3] == 1) & 
                                        (eccs[:,4] == 0)
                    ])
                    BC_count += len(eccs[(eccs[:,2] == 0) & 
                                        (eccs[:,3] == 1) & 
                                        (eccs[:,4] == 1)
                    ])
                    AC_count += len(eccs[(eccs[:,2] == 1) & 
                                        (eccs[:,3] == 0) & 
                                        (eccs[:,4] == 1)
                    ])
                    ABC_count += len(eccs[(eccs[:,2] == 1) & 
                                        (eccs[:,3] == 1) & 
                                        (eccs[:,4] == 1)
                    ])
                else:
                    A_count += len(eccs[(eccs[:,2] == 1) & 
                                        (eccs[:,3] == 0) 
                    ])
                    B_count += len(eccs[(eccs[:,2] == 0) & 
                                        (eccs[:,3] == 1) 
                    ])
                    AB_count += len(eccs[(eccs[:,2] == 1) & 
                                        (eccs[:,3] == 1) 
                    ])
        # print(biorep)
        # print([A_count, B_count, C_count,
        #       AB_count, BC_count, AC_count,
        #       ABC_count])
        list_of_lists_output.append([tolerance, biorep, A_count, B_count, C_count,
              AB_count, BC_count, AC_count,
              ABC_count])
    # print('now comparing all combined bioreps')
    bioreps_dict = {} ## now combine all techreps together and compare bioreps
    for biorep in all_bioreps:
        biorep_string = biorep[0][0:4]
        bioreps_dict[biorep_string] = []
        for i in range(56):
            samples_list = []
            for sample in biorep:
                if ecc_dict[sample][i].size != 0:
                    samples_list.append(ecc_dict[sample][i])
            if len(samples_list) != 0:
                samples_concatenated = np.concatenate((samples_list), axis=0) ## concatenate
            if samples_concatenated.size != 0:
                unique_samples_concatenated = np.unique(samples_concatenated, axis=0)
            else:
                unique_samples_concatenated = samples_concatenated
            bioreps_dict[biorep_string].append(unique_samples_concatenated)
    # print('getting data for all bioreps')
    overlap_bioreps = []
    for i in range(56):
        bioreps_list = []
        for biorep in bioreps_dict.keys():
            if bioreps_dict[biorep][i].size != 0:
                bioreps_list.append(bioreps_dict[biorep][i])
        bioreps_concatenated = np.concatenate((bioreps_list), axis=0)
        if bioreps_concatenated.size != 0:
            unique_bioreps_concatenated = np.unique(bioreps_concatenated, axis=0)
        else:
            unique_bioreps_concatenated = bioreps_concatenated
        if unique_bioreps_concatenated.size != 0:
            for biorep in bioreps_dict.keys():
                biorep_eccs = bioreps_dict[biorep][i]
                membership = []
                for ecc in unique_bioreps_concatenated:
                    if biorep_eccs.size != 0 and np.size(
                    biorep_eccs[np.logical_and(
                        np.isclose(biorep_eccs[:,0], ecc[0], atol=tolerance, rtol=0),
                        np.isclose(biorep_eccs[:,1], ecc[1], atol=tolerance, rtol=0)
                    )] ## add isclose or just == here, depending on the level of tolerance
                    ) != 0:
                        membership.append(1)
                    else:
                        membership.append(0)
                membership_column = np.reshape(np.array(membership), (-1, 1))
                unique_bioreps_concatenated = np.append(unique_bioreps_concatenated, membership_column, axis=1)
        overlap_bioreps.append(unique_bioreps_concatenated)
    # print('overlaps for all bioreps are...')
    # print('(order is A, B, C, AB, BC, AC, ABC)')
    count_1 = 0
    count_2 = 0
    count_3 = 0
    count_12 = 0
    count_23 = 0
    count_13 = 0
    count_123 = 0
    for i in range(56):
        eccs = overlap_bioreps[i]
        if eccs.size != 0:
            count_1 += len(eccs[(eccs[:,2] == 1) & 
                                (eccs[:,3] == 0) & 
                                (eccs[:,4] == 0)
            ])
            count_2 += len(eccs[(eccs[:,2] == 0) & 
                                (eccs[:,3] == 1) & 
                                (eccs[:,4] == 0)
            ])
            count_3 += len(eccs[(eccs[:,2] == 0) & 
                                (eccs[:,3] == 0) & 
                                (eccs[:,4] == 1)
            ])
            count_12 += len(eccs[(eccs[:,2] == 1) & 
                                (eccs[:,3] == 1) & 
                                (eccs[:,4] == 0)
            ])
            count_23 += len(eccs[(eccs[:,2] == 0) & 
                                (eccs[:,3] == 1) & 
                                (eccs[:,4] == 1)
            ])
            count_13 += len(eccs[(eccs[:,2] == 1) & 
                                (eccs[:,3] == 0) & 
                                (eccs[:,4] == 1)
            ])
            count_123 += len(eccs[(eccs[:,2] == 1) & 
                                (eccs[:,3] == 1) & 
                                (eccs[:,4] == 1)
            ])
        else:
            A_count += len(eccs[(eccs[:,2] == 1) & 
                                (eccs[:,3] == 0) 
            ])
            B_count += len(eccs[(eccs[:,2] == 0) & 
                                (eccs[:,3] == 1) 
            ])
            AB_count += len(eccs[(eccs[:,2] == 1) & 
                                (eccs[:,3] == 1) 
            ])
    # print('tech reps combined')
    # print([count_1, count_2, count_3,
    #       count_12, count_23, count_13,
    #       count_123])
    list_of_lists_output.append([tolerance, 'all', count_1, count_2, count_3,
          count_12, count_23, count_13,
          count_123])

with open(subdir+'/'+'venn_diagram_stats_out.tsv', 'w', newline = '') as out_file:
    w = csv.writer(out_file, delimiter = '\t')
    w.writerows(list_of_lists_output)