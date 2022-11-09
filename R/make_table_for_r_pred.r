library(tidyr)
library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

og_file = args[1]
gff_file = args[2]
te_file = args[3]
genome = args[4]
chromsizes_file = args[5]
gc_stats = args[6]
tmhmm_table = args[7]
signalp_table = args[8]
effectorp_table = args[9]
H3K27ac_per_og = args[10]
H3K27me3_per_og = args[11]
H3K36me3_per_og = args[12]
zhang_et_al_2019_complete_medium_expression_per_og = args[13]
zhang_et_al_2019_in_planta_expression_per_og = args[14]
eccdnacov_per_og = args[15]
methylation_per_og = args[16]
go_table = args[17]
pfam_table = args[18]
output_table = args[19]

## read in files
read_locations = function(file) {
  # list of elements to judge the distance of, in this case genes
  
  locations_file = file
  
  # list of elements to judge the distance to, in this case genes again
  
  
  # read in file
  # just genes and their location and orientation
  locations <- read.table(locations_file,header=FALSE,sep='\t',quote="")
  
  locations <- locations[locations$V3 == 'gene',]
  
  locations <- subset(locations, select = c(V1, V4, V5, V7, V9))
  
  colnames(locations) <- c('scaffold', 'start', 'end', 'orientation', 'id')
  
  locations$id <- gsub(".*ID=(.+);.*", "\\1", locations$id)
  
  locations$id <- gsub("^(.+);.*", "\\1", locations$id)
  
  locations$id <- gsub(".*gene:(.+);.*", "\\1", locations$id)
  
  locations$id <- gsub("^(.+);.*", "\\1", locations$id)
  
  locations$id <- paste(locations$id, 'T0', sep='')
  
  return(locations)
}

## read in comparison tes
read_comparison_te = function(file) {
  comparison_file = file
    
  comparison <- read.table(comparison_file, header=FALSE,sep='\t',quote="")
  
  comparison <- subset(comparison, select = c(V1, V4, V5, V6, V7, V9, V10))
  
  colnames(comparison) <- c('scaffold', 'start', 'end', 'identity', 'orientation', 'element', 'class')
  
  return(comparison)
}

## read in comparison genes
read_comparison_gene = function(file) {
  comparison_file = file
    
  comparison <- read.table(comparison_file, header=FALSE,sep='\t',quote="")

  comparison <- comparison[comparison$V3 == 'gene',]
  
  comparison <- subset(comparison, select = c(V1, V4, V5, V7, V9))
  
  colnames(comparison) <- c('scaffold', 'start', 'end', 'orientation', 'id')
  
  comparison$id <- gsub(".*ID=(.+);.*", "\\1", comparison$id)
  
  comparison$id <- gsub("^(.+);.*", "\\1", comparison$id)
  
  comparison$id <- gsub(".*gene:(.+);.*", "\\1", comparison$id)
  
  comparison$id <- gsub("^(.+);.*", "\\1", comparison$id)
  
  comparison$id <- paste(comparison$id, 'T0', sep='')
  
  return(comparison)
}

# add orthogroup information
add_ogs = function(genome_full_name, locations) {
  
  orthogroup_info <- data.frame(df_ogs[genome_full_name])
  
  rownames(orthogroup_info) <- rownames(df_ogs)
  
  orthogroup_info <- na.omit(orthogroup_info)
  
  colnames(orthogroup_info) <- c('protein')
  
  ## this is to deal with the separation later, need to know how big to make the data frame
  max_columns <- max(sapply(orthogroup_info$protein, strsplit_sapply))
  
  column_names <- as.character(seq(from = 1, to = max_columns, by = 1))
  
  orthogroup_info <- separate(data = orthogroup_info, col = protein, into = column_names, sep = ", ")
  
  og_list <- list()
  
  for (col in colnames(orthogroup_info)){
    og_list <- c(og_list, setNames(rownames(orthogroup_info[! is.na(orthogroup_info[[col]]), ]),
                                   na.omit(orthogroup_info[[col]])))
  }
  
  og_list_ordered <- og_list[locations$id]
  
  # replace null values with NAs
  og_list_ordered[sapply(og_list_ordered, is.null)] <- NA
  
  return(as.vector(unlist(og_list_ordered)))
}

strsplit_sapply = function(gene_string){
  return(length(strsplit(gene_string, ', ')[[1]]))
}


read_chromsizes = function(file) {
  chromsizes_file = file
  
  chromsizes <- read.table(chromsizes_file,header=FALSE,sep='\t',quote="")
  
  colnames(chromsizes) <- c('scaffold','size')
  
  return(chromsizes)
}

## check if there is enough space around the gene to call
enough_space = function(x, distance) {
  scaffold = x['scaffold']
  start = as.numeric(x['start'])
  end = as.numeric(x['end'])
  scaffold_size = chromsizes$size[chromsizes$scaffold==scaffold]
  if (start >= distance | end <= (scaffold_size - distance)){ ## check if there are enough sequenced base pairs to call a FALSE on at least one side
      return(TRUE)
  } else {
    return(FALSE)
  }
}

## check if a te is nearby or not
nearby_te_per_class = function(x, class, distance) {
  scaffold = x['scaffold']
  start = as.numeric(x['start'])
  end = as.numeric(x['end'])
  enough_space_te = x['enough_space_te']
  nearby_te = subset_per_class_per_scaffold[subset_per_class_per_scaffold$end >= start - distance
                      & subset_per_class_per_scaffold$start <= end + distance,]
  if (dim(nearby_te)[1] == 0){ # if empty
    if (enough_space_te == TRUE){ ## check if there are enough sequenced base pairs to call a FALSE
      return(FALSE)
    } else {
      return(NA) ## otherwise we dont know
    }
  } else{
    return(TRUE)
  }
}

nearby_gene = function(x, class, distance) {
  scaffold = x['scaffold']
  start = as.numeric(x['start'])
  end = as.numeric(x['end'])
  enough_space_gene = x['enough_space_gene']
  id = x['id']
  nearby_genes = subset_genes_per_scaffold[subset_genes_per_scaffold$end >= start - distance
                      & subset_genes_per_scaffold$start <= end + distance
                      & subset_genes_per_scaffold$id != id,] ## make sure you dont match to the same gene each time
  if (dim(nearby_genes)[1] == 0){ # if empty
    if (enough_space_gene == TRUE){ ## check if there are enough sequenced base pairs to call a FALSE
      return(FALSE)
    } else {
      return(NA) ## otherwise we dont know
    }
  } else{
    return(TRUE)
  }
}

## read in ogs
df_ogs <- data.frame(fread(og_file))

rownames(df_ogs) <- df_ogs$Orthogroup

df_ogs <- subset(df_ogs, select = -c(Orthogroup))

df_gene_info <- data.frame(matrix(ncol = 11, nrow = 0))

colnames(df_gene_info) <- c("scaffold", "start", "end", "orientation", "id", "orthogroups", 
                            "enough_space_te", "enough_space_gene",
                              'any_te',
                            "gene_nearby",
                            "genome"
                            )

locations <- read_locations(gff_file)
comparison_tes <- read_comparison_te(te_file)
comparison_genes <- read_comparison_gene(gff_file)

# add og information to locations AND comparison
locations$orthogroups <- suppressWarnings(add_ogs(genome, locations))
# drop genes without ogs
locations <- locations[! (is.na(locations$orthogroups)),]
comparison_genes$orthogroups <- suppressWarnings(add_ogs(genome, comparison_genes))

chromsizes <- read_chromsizes(chromsizes_file)

distance_te = 5000

locations$enough_space_te = apply(locations,1,enough_space,distance=distance_te)

distance_gene = 1000

locations$enough_space_gene = apply(locations,1,enough_space,distance=distance_gene)

## make sure locations is sorted by scaffold

locations = locations[order(locations$scaffold),]

class = "any_te"
if (class == "any_te"){
  subset_per_class = comparison_tes
} else {
  subset_per_class = comparison_tes[comparison_tes$class == class,]
}
if (dim(subset_per_class)[1] != 0){
  class_column = c()
  for (scaffold in unique(locations$scaffold)){
    locations_per_scaffold = locations[locations$scaffold == scaffold,]
    subset_per_class_per_scaffold = subset_per_class[subset_per_class$scaffold == scaffold,]
    if (dim(subset_per_class_per_scaffold)[1] != 0){
      locations_per_scaffold[[class]] = apply(locations_per_scaffold,1,nearby_te_per_class,distance=distance_te)
    } else {
      locations_per_scaffold[[class]][locations_per_scaffold$enough_space_te == TRUE] = FALSE
    }
    class_column = c(class_column, locations_per_scaffold[[class]])
  }
  locations[[class]] = class_column
} else {
  locations[[class]][locations$enough_space_te == TRUE] = FALSE ## if there are no TEs of that class in that genome then set all to FALSE
}


gene_column = c()

for (scaffold in unique(locations$scaffold)){
  locations_per_scaffold = locations[locations$scaffold == scaffold,]
  subset_genes_per_scaffold = comparison_genes[comparison_genes$scaffold == scaffold,]
  # there will always be at least one gene...
  gene_column = c(gene_column, apply(locations_per_scaffold,1,nearby_gene,distance=distance_gene))
}

locations['gene_nearby'] = gene_column

locations$genome <- genome

locations <- locations[,colnames(df_gene_info)]

df_gene_info <- rbind(df_gene_info, locations)

df_gc_stats <- data.frame(fread(gc_stats,sep='\t'))

colnames(df_gc_stats) <- c("id", "gene_gc", "flanking_1kb_gc")

df_gc_stats <- subset(df_gc_stats, select =  c("id", "gene_gc", "flanking_1kb_gc"))

df_gc_stats$id <- gsub(".*ID=(.+);.*", "\\1", df_gc_stats$id)

df_gc_stats$id <- gsub("^(.+);.*", "\\1", df_gc_stats$id)

df_gc_stats$id <- gsub(".*gene:(.+);.*", "\\1", df_gc_stats$id)

df_gc_stats$id <- gsub(".*gene:(.+).*", "\\1", df_gc_stats$id)

df_gc_stats$id <- paste(df_gc_stats$id, 'T0', sep='')

df_gene_info <- merge(df_gene_info, df_gc_stats, by = "id")

rownames(df_gene_info) <- df_gene_info$id

## add gene lengths

df_gene_info$lengths <- df_gene_info$end - df_gene_info$start

df_tm <- data.frame(fread(tmhmm_table), na.string='')

colnames(df_tm) <- c("id", "tm")

df_tm <- subset(df_tm, select = c("id", "tm"))

df_gene_info <- merge(df_gene_info, df_tm, by = "id")

rownames(df_gene_info) <- df_gene_info$id



df_signalp <- data.frame(fread(signalp_table), na.string='')

colnames(df_signalp) <- c("id", "signalp")

df_signalp <- subset(df_signalp, select = c("id", "signalp"))

df_gene_info <- merge(df_gene_info, df_signalp, by = "id")

rownames(df_gene_info) <- df_gene_info$id


df_effectorp <- data.frame(fread(effectorp_table), na.string='')

colnames(df_effectorp) <- c("id", "effectorp")

df_effectorp <- subset(df_effectorp, select = c("id", "effectorp"))

df_gene_info <- merge(df_gene_info, df_effectorp, by = "id", all.x = TRUE)

rownames(df_gene_info) <- df_gene_info$id

df_gene_info$effectorp[is.na(df_gene_info$effectorp)] <- FALSE

signal_files <- c(
  H3K27ac_per_og,
  H3K27me3_per_og,
  H3K36me3_per_og,
  zhang_et_al_2019_complete_medium_expression_per_og,
  zhang_et_al_2019_in_planta_expression_per_og,
  eccdnacov_per_og,
  methylation_per_og
  )

signal_names <- c('H3K27ac',
    'H3K27me3',
    'H3K36me3',
    'cm_expression',
    'ip_expression',
    'eccdna_cov',
    'methylation')

for (i in seq(length(signal_names))) {
  signal_file = signal_files[i]
  signal_name = signal_names[i]
  print(signal_file)
  print(signal_name)
  
  df_signal <- data.frame(fread(signal_file), na.string='')

  colnames(df_signal) <- c("orthogroups", signal_name)
  
  df_signal <- subset(df_signal, select = c("orthogroups", signal_name))
  
  df_gene_info <- merge(df_gene_info, df_signal, by = "orthogroups")
  
  rownames(df_gene_info) <- df_gene_info$id

}

df_go <- data.frame(fread(go_table), na.string='')

colnames(df_go) <- c("id", "go")

df_go <- subset(df_go, select = c("id", "go"))

df_gene_info <- merge(df_gene_info, df_go, by = "id", all.x = TRUE)

rownames(df_gene_info) <- df_gene_info$id

df_gene_info$go[is.na(df_gene_info$go)] <- FALSE

df_pfam <- data.frame(fread(pfam_table), na.string='')

colnames(df_pfam) <- c("id", "pfam")

df_pfam <- subset(df_pfam, select = c("id", "pfam"))

df_gene_info <- merge(df_gene_info, df_pfam, by = "id", all.x = TRUE)

rownames(df_gene_info) <- df_gene_info$id

df_gene_info$pfam[is.na(df_gene_info$pfam)] <- FALSE

## and drop all NAs bc random forest doesnt work with them
df_gene_info <- na.omit(df_gene_info)

output.file <- file(output_table, "wb") ## to get unix line endings

write.table(df_gene_info, file = output.file, quote = FALSE, col.names = TRUE, row.names = TRUE, sep=',')

close(output.file)