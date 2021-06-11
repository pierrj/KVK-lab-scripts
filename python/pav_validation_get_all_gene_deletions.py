import csv
from os import listdir
import sys
import os
import shutil


directory_orthogroups = str(sys.argv[1])
fungap_out = str(sys.argv[2])
directory_gffs = str(sys.argv[3])
directory_sizes = str(sys.argv[4])
output_directory = str(sys.argv[5])

list_of_files = sorted(listdir(directory_orthogroups))
list_of_files.remove('guy11_fungap_out_prot__v__Neurospora_crassa.NC12.pep.all.tsv')
list_of_files.remove('guy11_fungap_out_prot__v__GCA_002368485.1_ASM236848v1_fungap_out_prot.tsv')

guy11_gene_list = []
with open(fungap_out, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        try:
            if row[2] == 'gene':
                guy11_gene_list.append(row[8][3:13]+'.t1')
        except IndexError:
            pass

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

## skip blastn step
genome_dict_actually_missing = genomes_dict

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
                guy11_gene_list_long.append([row[0], row[3], row[4], row[8][3:13]+'.t1'])
                guy11_gene_list_short.append(row[8][3:13]+'.t1')
                #gene_list.append([row[0], row[3], row[4], row[8][2:11]])
        except IndexError:
            pass
guy11_neighbors = {}
for i in range(len(guy11_gene_list_long)):
    if i-1 < 0 or guy11_gene_list_long[i-1][0] != guy11_gene_list_long[i][0]:
        gene_1 = 'none'
    else:
        gene_1 = guy11_gene_list_long[i-1][3]
    gene_2 = guy11_gene_list_long[i][3]
    if i == len(guy11_gene_list_long)-1 or guy11_gene_list_long[i+1][0] != guy11_gene_list_long[i][0]:
        gene_3 = 'none'
    else:
        gene_3 = guy11_gene_list_long[i+1][3]
    guy11_neighbors[gene_2] = [gene_1,gene_3]

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
                    genome_long_list.append([row[0], row[3], row[4], row[8][3:13]+'.t1'])
                    genome_short_list.append(row[8][3:13]+'.t1')
                    #gene_list.append([row[0], row[3], row[4], row[8][2:11]])
            except IndexError:
                pass
    genome_neighbors_dict['no_neighbors'] = []
    for i in range(len(genome_long_list)):
        gene_2 = genome_long_list[i][3]
        if i-1 < 0 or genome_long_list[i-1][0] != genome_long_list[i][0]:
            genome_neighbors_dict['no_neighbors'].append(gene_2)
            continue ## maybe just exclude these
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

for i in range(len(list_of_genomes)):
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
                    #gene_list.append([row[0], row[3], row[4], row[8][2:11]])
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
translation_dictionary_guy11_key = {}
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
            if len(gene_guy11.split(', ')) > 1:
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
            translation_dictionary_pergenome['genome_genes_without_guy11_orthologs'].append(genome_gene)
    for k in range(len(guy11_gene_list_short)):
        guy11_gene = guy11_gene_list_short[k]
        if guy11_gene not in found_guy11_genes:
            translation_dictionary_pergenome['guy11_genes_without_genome_orthologs'].append(guy11_gene)

lost_genes = []
for i in range(len(guy11_gene_list_short)):
    gene = guy11_gene_list_short[i]
    if pav_final[gene] != 162:
        lost_genes.append(gene)

gene_deletion_dict = {}
deleted_gene_locs = {}
for l in range(len(lost_genes)):
    missing_gene_of_interest = lost_genes[l]
    missing_gene_of_interest_loc = 0
    gene_deletion_dict[missing_gene_of_interest] = []
    for k in range(len(guy11_gene_list_long)):
        if guy11_gene_list_long[k][3] == missing_gene_of_interest:
            missing_gene_of_interest_loc = guy11_gene_list_long[k]
    deleted_gene_locs[missing_gene_of_interest] = missing_gene_of_interest_loc
    genomes_missing_gene = []
    for genome in genome_dict_actually_missing.keys():
        for i in range(len(genome_dict_actually_missing[genome])):
            if genome_dict_actually_missing[genome][i] == missing_gene_of_interest:
                genomes_missing_gene.append(genome)
    neighbors = guy11_neighbors[missing_gene_of_interest]
    neighbors_translated = []
    for i in range(len(genomes_missing_gene)):
        genome = genomes_missing_gene[i]
        genome_translation_dict = translation_dictionary_guy11_key[genome]
        neighbor_5_translated = 0
        neighbor_3_translated = 0
        try:
            neighbor_5_translated = genome_translation_dict[neighbors[0]]
        except KeyError:
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
        except KeyError:
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
    neighbors_translated_with_locs = []
    for i in range(len(neighbors_translated)):
        genome = neighbors_translated[i][0]
        neighbor_5_translated = neighbors_translated[i][1]
        neighbor_3_translated = neighbors_translated[i][2]
        gene_list_long = gene_list_long_dict[genome]
        neighbor_5_translated_with_loc = 0
        neighbor_3_translated_with_loc = 0
        for k in range(len(gene_list_long)):
            if gene_list_long[k][3] == neighbor_5_translated:
                neighbor_5_translated_with_loc = gene_list_long[k]
            if gene_list_long[k][3] == neighbor_3_translated:
                neighbor_3_translated_with_loc = gene_list_long[k]
        neighbors_translated_with_locs.append([genome, neighbor_5_translated_with_loc, neighbor_3_translated_with_loc])
    for i in range(len(neighbors_translated_with_locs)):
        if (0 not in neighbors_translated_with_locs[i] and
           abs(int(neighbors_translated_with_locs[i][1][3][5:10]) - int(neighbors_translated_with_locs[i][2][3][5:10])) < 2 and
           neighbors_translated_with_locs[i][1][0] == neighbors_translated_with_locs[i][2][0]):
            gene_deletion_dict[missing_gene_of_interest].append(neighbors_translated_with_locs[i])

list_of_genome_sizes = sorted(listdir(directory_sizes))
genomesize_dict = {}
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
os.mkdir(output_directory)
reference_name = 'guy11'
for gene in gene_deletion_dict.keys():
    if len(gene_deletion_dict[gene]) > 0:
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
            start_ref = int(gene_loc[1])
            end_ref = int(gene_loc[2])
            start_quer = min([int(deletion[1][1]), int(deletion[1][2]), int(deletion[2][1]), int(deletion[2][2])])
            end_quer = max([int(deletion[1][1]), int(deletion[1][2]), int(deletion[2][1]), int(deletion[2][2])])
            if start_ref - 20000 > 0:
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
                w = csv.writer(output_csv, delimiter = '\t') ## lineterminator="\n" might fix issues with windows
                w.writerow([gene_loc[0], start_ref, end_ref])
            with open(output_directory+'/'+gene+'/'+outdir_for_comparison + "/" + query_name + '_0.bed', 'w', newline = '') as output_csv:
                w = csv.writer(output_csv, delimiter = '\t') ## lineterminator="\n" might fix issues with windows
                w.writerow([deletion[1][0], start_quer, end_quer])
            with open(output_directory+'/'+gene+'/'+outdir_for_comparison+"/"+reference_name+'_geneloc.bed', 'w', newline = '') as output_csv:
                w = csv.writer(output_csv, delimiter = '\t')
                w.writerow([gene_loc[0], int(gene_loc[1])-start_ref, int(gene_loc[2])-start_ref])