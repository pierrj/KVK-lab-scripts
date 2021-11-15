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
from os import listdir
import sys
import os
import shutil

## USAGE ##
# this script takes input from orthofinder (i.e. guy11 vs all other genomes) and outputs a list of genes that are deleted from the location and the location of where that gene deletion is
# location of gene deletion is based off the presence of neighboring gene equivalents in other genomes (according to orthogroups from orthofinder)
# "directory_orthogroups" - directory of orthogroup comparison results (i.e. guy11 vs all other genomes)
# "fungap_out" - gene gff file for reference
# "directory_gffs" - directory containing gff files for all other genomes
# "directory_sizes" - directory of genome size files (scaffold names in one column, scaffold length in the other column) for all other genomes
# "output_directory" - directory for outputs of bed files where putative deletion locations are

directory_orthogroups = str(sys.argv[1])
fungap_out = str(sys.argv[2])
directory_gffs = str(sys.argv[3])
directory_sizes = str(sys.argv[4])
output_directory = str(sys.argv[5])

# remove ncrassa from output as well as redudant guy11 genome
list_of_files = sorted(listdir(directory_orthogroups))
list_of_files.remove('guy11_fungap_out_prot__v__Neurospora_crassa.NC12.pep.all.tsv')
list_of_files.remove('guy11_fungap_out_prot__v__GCA_002368485.1_ASM236848v1_fungap_out_prot.tsv')

# read in list of genes in guy11
guy11_gene_list = []
with open(fungap_out, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        try:
            if row[2] == 'gene':
                guy11_gene_list.append(row[8][3:13]+'.t1')
        except IndexError:
            pass

# read in orthofinder output into dictionary 
genomes_dict = {}
for i in range(len(list_of_files)):
    if list_of_files[i] == 'guy11_fungap_out_prot__v__Neurospora_crassa.NC12.pep.all.tsv':
        continue
    genome = list_of_files[i][26:-20]
    found_genes_list = []
    with open(directory_orthogroups + '/' + list_of_files[i], newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        next(file_reader)
        for row in file_reader:
            gene = row[1]
            if len(gene.split(', ')) > 1:
                for ii in range(len(gene.split(', '))):
                    split_gene = gene.split(', ')[ii]
                    found_genes_list.append(split_gene)
            else:
                found_genes_list.append(gene)
    genomes_dict[genome] = list(set(guy11_gene_list) - set(found_genes_list))

## skip blastn validation step from other pav analyses
genome_dict_actually_missing = genomes_dict

## count how many times each orthogroup is missing from other genomes based off orthofinder output
pav_final = {}
for i in range(len(guy11_gene_list)):
    gene = guy11_gene_list[i]
    pav_final[gene] = 162
for genome in genome_dict_actually_missing.keys():
    for i in range(len(genome_dict_actually_missing[genome])):
        pav_final[genome_dict_actually_missing[genome][i]] -= 1

gff_directory = directory_gffs 
list_of_genomes = sorted(listdir(gff_directory))
list_of_genomes.remove('guy11_fungap_out.gff3')
list_of_genomes.remove('GCA_002368485.1_ASM236848v1_fungap_out.gff3')
orthologs_directory = directory_orthogroups
list_of_orthologs = sorted(listdir(orthologs_directory))
list_of_orthologs.remove('guy11_fungap_out_prot__v__Neurospora_crassa.NC12.pep.all.tsv')
list_of_orthologs.remove('guy11_fungap_out_prot__v__GCA_002368485.1_ASM236848v1_fungap_out_prot.tsv')

guy11_gene_list_long = []
guy11_gene_list_short = []
with open(gff_directory + '/' + 'guy11_fungap_out.gff3', newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        try:
            if row[2] == 'gene':
                guy11_gene_list_long.append([row[0], row[3], row[4], row[8][3:13]+'.t1']) # this list contains the name and location of genes
                guy11_gene_list_short.append(row[8][3:13]+'.t1') # this list just contains the name of genes
        except IndexError: # some issues with reading in files sometimes...
            pass
guy11_neighbors = {} # for each gene in guy11, find which genes are direct neighbors of this gene, keys are the genes
for i in range(len(guy11_gene_list_long)):
    if i-1 < 0 or guy11_gene_list_long[i-1][0] != guy11_gene_list_long[i][0]: # upstream gene, making sure that it isnt the first gene, and that the upstream gene is on hte same scaffold
        gene_1 = 'none'
    else:
        gene_1 = guy11_gene_list_long[i-1][3]
    gene_2 = guy11_gene_list_long[i][3]
    if i == len(guy11_gene_list_long)-1 or guy11_gene_list_long[i+1][0] != guy11_gene_list_long[i][0]: # same as before but with downstream
        gene_3 = 'none'
    else:
        gene_3 = guy11_gene_list_long[i+1][3]
    guy11_neighbors[gene_2] = [gene_1,gene_3]

# get all three as before for guy11, but this time for each other genome
gene_list_long_dict = {}
gene_list_short_dict = {}
gene_neighbors_dict = {}

for i in range(len(list_of_genomes)):
    gff_file = list_of_genomes[i]
    genome = gff_file[:-16]
    gene_list_long_dict[genome] = []
    gene_list_short_dict[genome] = []
    gene_neighbors_dict[genome] = {}
    genome_long_list = gene_list_long_dict[genome]
    genome_short_list = gene_list_short_dict[genome]
    genome_neighbors_dict = gene_neighbors_dict[genome]
    with open(gff_directory + '/' + gff_file, newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        for row in file_reader:
            try:
                if row[2] == 'gene':
                    genome_long_list.append([row[0], row[3], row[4], row[8][3:13]+'.t1']) # read in gene locs
                    genome_short_list.append(row[8][3:13]+'.t1') # read in gene names
            except IndexError:
                pass
    genome_neighbors_dict['no_neighbors'] = [] # make list of genes without neighbors
    for i in range(len(genome_long_list)): ## store gene neighbors as before
        gene_2 = genome_long_list[i][3]
        if i-1 < 0 or genome_long_list[i-1][0] != genome_long_list[i][0]:
            genome_neighbors_dict['no_neighbors'].append(gene_2)
            continue
        else:
            gene_1 = genome_long_list[i-1][3]
        if i == len(genome_long_list)-1 or genome_long_list[i+1][0] != genome_long_list[i][0]:
            genome_neighbors_dict['no_neighbors'].append(gene_2)
            continue
        else:
            gene_3 = genome_long_list[i+1][3]
        genome_neighbors_dict[gene_2] = [gene_1,gene_3]

gene_list_long_nones_dict = {}
gene_list_short_nones_dict = {}
gene_neighbors_nones_dict = {}

for i in range(len(list_of_genomes)): ## deal with genes taht are missing neighbros and write them to gene_neighbors_dict
    gff_file = list_of_genomes[i]
    genome = gff_file[:-16]
    gene_list_long_nones_dict[genome] = []
    gene_list_short_nones_dict[genome] = []
    gene_neighbors_nones_dict[genome] = {}
    genome_long_list = gene_list_long_nones_dict[genome]
    genome_short_list = gene_list_short_nones_dict[genome]
    genome_neighbors_dict = gene_neighbors_nones_dict[genome]
    with open(gff_directory + '/' + gff_file, newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        for row in file_reader:
            try:
                if row[2] == 'gene':
                    genome_long_list.append([row[0], row[3], row[4], row[8][3:13]+'.t1'])
                    genome_short_list.append(row[8][3:13]+'.t1')
            except IndexError:
                pass
    for i in range(len(genome_long_list)):
        gene_2 = genome_long_list[i][3]
        if i-1 < 0 or genome_long_list[i-1][0] != genome_long_list[i][0]:
            gene_1 = 'none'
        else:
            gene_1 = genome_long_list[i-1][3]
        if i == len(genome_long_list)-1 or genome_long_list[i+1][0] != genome_long_list[i][0]:
            gene_3 = 'none'
        else:
            gene_3 = genome_long_list[i+1][3]
        genome_neighbors_dict[gene_2] = [gene_1,gene_3]

## this one takes a while
translation_dictionary_guy11_key = {} ## make a big dictionary that uses the orthofinder orthogroup information to make an easy way to translate between guy11 genes and other genome genes
# guy11 is the key here so it lets you translate guy11 genes to other genome genes
for i in range(len(list_of_orthologs)):
    ortholog_file = list_of_orthologs[i]
    genome = ortholog_file[26:-20]
    translation_dictionary_guy11_key[genome] = {}
    translation_dictionary_pergenome = translation_dictionary_guy11_key[genome]
    found_genome_genes = []
    found_guy11_genes = []
    translation_dictionary_pergenome['guy11_not_one_to_one'] = []
    translation_dictionary_pergenome['genome_not_one_to_one'] = []
    translation_dictionary_pergenome['genome_genes_without_guy11_orthologs'] = []
    translation_dictionary_pergenome['guy11_genes_without_genome_orthologs'] = []
    with open(orthologs_directory + '/'+ ortholog_file, newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        next(file_reader)
        for row in file_reader:
            gene_guy11 = row[1]
            gene_genome = row[2]
            if len(gene_guy11.split(', ')) > 1: ## deal with not one to one orthogroups which are excluded from the analysis
                for ii in range(len(gene_guy11.split(', '))):
                    split_gene = gene_guy11.split(', ')[ii]
                    translation_dictionary_pergenome['guy11_not_one_to_one'].append(split_gene)
                    found_guy11_genes.append(split_gene)
            if len(gene_genome.split(', ')) > 1:
                for ii in range(len(gene_genome.split(', '))):
                    split_gene = gene_genome.split(', ')[ii]
                    translation_dictionary_pergenome['genome_not_one_to_one'].append(split_gene)
                    found_genome_genes.append(split_gene)
            else:
                translation_dictionary_pergenome[gene_guy11] = gene_genome
                found_guy11_genes.append(gene_guy11)
                found_genome_genes.append(gene_genome)
    for k in range(len(gene_list_short_dict[genome])):
        genome_gene = gene_list_short_dict[genome][k]
        if genome_gene not in found_genome_genes:
            translation_dictionary_pergenome['genome_genes_without_guy11_orthologs'].append(genome_gene) ## deal with genes without guy11 orthologs
    for k in range(len(guy11_gene_list_short)):
        guy11_gene = guy11_gene_list_short[k]
        if guy11_gene not in found_guy11_genes:
            translation_dictionary_pergenome['guy11_genes_without_genome_orthologs'].append(guy11_gene) ## deal with guy11 genes without orthologs in the genome

## list of genes lost in genome based off pav score
lost_genes = []
for i in range(len(guy11_gene_list_short)):
    gene = guy11_gene_list_short[i]
    if pav_final[gene] != 162:
        lost_genes.append(gene)


# write deleted genes of interest and their expected locations in the genome
gene_deletion_dict = {}
deleted_gene_locs = {}
for l in range(len(lost_genes)): ## loop through lost genes
    missing_gene_of_interest = lost_genes[l]
    missing_gene_of_interest_loc = 0
    gene_deletion_dict[missing_gene_of_interest] = []
    for k in range(len(guy11_gene_list_long)):
        if guy11_gene_list_long[k][3] == missing_gene_of_interest:
            missing_gene_of_interest_loc = guy11_gene_list_long[k]
    deleted_gene_locs[missing_gene_of_interest] = missing_gene_of_interest_loc # find location of lost gene
    genomes_missing_gene = [] # add list of genomes that are missing the gene
    for genome in genome_dict_actually_missing.keys():
        for i in range(len(genome_dict_actually_missing[genome])):
            if genome_dict_actually_missing[genome][i] == missing_gene_of_interest:
                genomes_missing_gene.append(genome)
    neighbors = guy11_neighbors[missing_gene_of_interest] # get the neighbors of the missing genes
    neighbors_translated = []
    for i in range(len(genomes_missing_gene)): # translate the neighbors to what the gene is in all other genomes
        genome = genomes_missing_gene[i]
        genome_translation_dict = translation_dictionary_guy11_key[genome]
        neighbor_5_translated = 0
        neighbor_3_translated = 0
        try:
            neighbor_5_translated = genome_translation_dict[neighbors[0]]
        except KeyError: # if not 1:1 ortholog or no neighbor, drop the gene of interest
            if neighbors[0] in genome_translation_dict['guy11_genes_without_genome_orthologs']:
                neighbor_5_translated = 'no_genome_equivalent'
            if neighbors[0] in genome_translation_dict['genome_not_one_to_one']:
                neighbor_5_translated = 'one_of_multiple_orthologs_in_genome'
            if neighbors[0] in genome_translation_dict['guy11_not_one_to_one']:
                neighbor_5_translated = 'multiple_orthologs_in_guy11'
            if neighbors[0] == 'none':
                neighbor_5_translated = 'neighbor_missing'
            if neighbor_5_translated == 0:
                print('not sure what happened')
        try:
            neighbor_3_translated = genome_translation_dict[neighbors[1]]
        except KeyError: # same thing for second neighbor
            if neighbors[1] in genome_translation_dict['guy11_genes_without_genome_orthologs']:
                neighbor_3_translated = 'no_genome_equivalent'
            if neighbors[1] in genome_translation_dict['genome_not_one_to_one']:
                neighbor_3_translated = 'one_of_multiple_orthologs_in_genome'
            if neighbors[1] in genome_translation_dict['guy11_not_one_to_one']:
                neighbor_3_translated = 'multiple_orthologs_in_guy11'
            if neighbors[1] == 'none':
                neighbor_3_translated = 'neighbor_missing'
            if neighbor_3_translated == 0:
                print('not sure what happened')
        neighbors_translated.append([genome, neighbor_5_translated, neighbor_3_translated])
    neighbors_translated_with_locs = [] ## get the location of the neighbors in the genome
    for i in range(len(neighbors_translated)):
        genome = neighbors_translated[i][0]
        neighbor_5_translated = neighbors_translated[i][1] # first and second neighbor
        neighbor_3_translated = neighbors_translated[i][2]
        gene_list_long = gene_list_long_dict[genome]
        neighbor_5_translated_with_loc = 0
        neighbor_3_translated_with_loc = 0
        for k in range(len(gene_list_long)): # grab location of neighbors in genome
            if gene_list_long[k][3] == neighbor_5_translated:
                neighbor_5_translated_with_loc = gene_list_long[k]
            if gene_list_long[k][3] == neighbor_3_translated:
                neighbor_3_translated_with_loc = gene_list_long[k]
        neighbors_translated_with_locs.append([genome, neighbor_5_translated_with_loc, neighbor_3_translated_with_loc])
    for i in range(len(neighbors_translated_with_locs)): # make sure that neighbors in guy11 and each others neighbors in the genome, basically just a single gene location
        if (0 not in neighbors_translated_with_locs[i] and
           abs(int(neighbors_translated_with_locs[i][1][3][5:10]) - int(neighbors_translated_with_locs[i][2][3][5:10])) < 2 and ## genes are neighbors if their gene names are right next to each other
           neighbors_translated_with_locs[i][1][0] == neighbors_translated_with_locs[i][2][0]): ## make sure they aren the same gene
            gene_deletion_dict[missing_gene_of_interest].append(neighbors_translated_with_locs[i])

list_of_genome_sizes = sorted(listdir(directory_sizes)) ## needs to be sorted to work on linux...
genomesize_dict = {} ## get genome sizes files for all genomes
for i in range(len(list_of_genome_sizes)):
    size = list_of_genome_sizes[i]
    genome = size[:-11]
    genomesize_dict[genome] = {}
    with open(directory_sizes+'/'+size, newline = '') as file:
            file_reader = csv.reader(file, delimiter = '\t')
            for row in file_reader:
                genomesize_dict[genome][row[0]] = int(row[1])
genomesize_dict['guy11'] = {}
with open(directory_sizes+'/guy11_genome_baoetal2017.fasta.genomesize', newline = '') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        for row in file_reader:
            genomesize_dict['guy11'][row[0]] = int(row[1])

if os.path.isdir(output_directory):
    shutil.rmtree(output_directory)
os.mkdir(output_directory) # make output directory
reference_name = 'guy11'
for gene in gene_deletion_dict.keys():
    if len(gene_deletion_dict[gene]) > 0: # make an output directory for each guy11 deleted gene
        if os.path.isdir(output_directory+'/'+gene):
            shutil.rmtree(output_directory+'/'+gene)
        os.mkdir(output_directory+'/'+gene)
        gene_loc = deleted_gene_locs[gene]
        reference_genomesize = genomesize_dict[reference_name]
        for i in range(len(gene_deletion_dict[gene])):
            deletion = gene_deletion_dict[gene][i]
            query_name = deletion[0]
            outdir_for_comparison = reference_name+'_v_'+deletion[0]
            os.mkdir(output_directory+'/'+gene+'/'+outdir_for_comparison)
            query_genomesize = genomesize_dict[query_name]
            start_ref = int(gene_loc[1]) # get location of deleted gene
            end_ref = int(gene_loc[2])
            start_quer = min([int(deletion[1][1]), int(deletion[1][2]), int(deletion[2][1]), int(deletion[2][2])]) # get location of putative deletion
            end_quer = max([int(deletion[1][1]), int(deletion[1][2]), int(deletion[2][1]), int(deletion[2][2])])
            if start_ref - 20000 > 0: # extend outward 20k base pairs on either side for alignment, dealing with scaffold ends appropriately
                start_ref = start_ref - 20000
            else:
                start_ref = 1
            if end_ref + 20000 < reference_genomesize[gene_loc[0]]:
                end_ref = end_ref + 20000
            else:
                end_ref = reference_genomesize[gene_loc[0]]
            if start_quer - 20000 > 0:
                start_quer = start_quer - 20000
            else:
                start_quer = 1
            if end_quer + 20000 < query_genomesize[deletion[1][0]]:
                end_quer = end_quer + 20000
            else:
                end_quer = query_genomesize[deletion[1][0]]
            with open(output_directory+'/'+gene+'/'+outdir_for_comparison + "/" + reference_name + '_0.bed', 'w', newline = '') as output_csv:
                w = csv.writer(output_csv, delimiter = '\t') 
                w.writerow([gene_loc[0], start_ref, end_ref]) # write bed file of deleted gene location in guy11 + 20k bp flanking
            with open(output_directory+'/'+gene+'/'+outdir_for_comparison + "/" + query_name + '_0.bed', 'w', newline = '') as output_csv:
                w = csv.writer(output_csv, delimiter = '\t')
                w.writerow([deletion[1][0], start_quer, end_quer]) # write bed file of putative location of deleted gene in query genome
            with open(output_directory+'/'+gene+'/'+outdir_for_comparison+"/"+reference_name+'_geneloc.bed', 'w', newline = '') as output_csv:
                w = csv.writer(output_csv, delimiter = '\t')
                w.writerow([gene_loc[0], int(gene_loc[1])-start_ref, int(gene_loc[2])-start_ref]) # write location of where the gene should be in the alignment between the two previous files