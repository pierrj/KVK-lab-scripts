library(data.table)


## read in dfs from random forest gene info tables for wheat and rice blast

df_gene_info_wheat_blast <-data.frame(fread('gene_info.cross_host.wheat_blast.txt', na.strings = ''))

df_gene_info_rice_blast <- data.frame(fread('gene_info.full_model.rice_blast.txt', na.strings = ''))


twofive_quant = function(vector){
  return(unname(quantile(vector, 0.25)))
}

sevenfive_quant = function(vector){
  return(unname(quantile(vector,0.75)))
}

# get a bunch of stats for an input vector
get_stats = function(vector, stats) {
  
  output_stats <- c()
  
  for (stat in stats){
    
    output_number <- stat(vector)
    
    output_stats <- c(output_stats, output_number)
  }
  
  return(output_stats)
}

permutation_test = function(vector1, vector2, stat, replicates){
    observed = stat(vector1) - stat(vector2)
    combined = c(vector1, vector2)

    null_dsn <- replicate(replicates,
                          {
                            sample_one_index <- sample(length(combined), length(vector1), replace = FALSE)
                            sample_two_index <- seq(length(combined))[!seq(length(combined))%in%sample_one_index]

                            sample_one <- combined[sample_one_index]
                            sample_two <- combined[sample_two_index]
                            stat(sample_one) - stat(sample_two)
                          })

    ## this is the two tailed p-value
    p_value <- mean(abs(observed)<=abs(null_dsn))

    if (p_value == 0){
      return(paste("<", as.character(1/replicates), sep=''))
    } else {
      return(as.character(p_value))
    }
}

# permutation_test = function(vector1, vector2, stat, replicates){
#   return('NA')
# }

get_all_pvalues = function(vector1, vector2, stats, replicates) {
  p_values <- c()
  
  for (stat in stats) {
    p_values <- c(p_values, permutation_test(vector1, vector2, stat, replicates))
  }
  
  return(p_values)
}

## params for which statistics to generate

stats <- 
  c(mean, median, sd, twofive_quant, sevenfive_quant)

stats_labels <- 
  c("mean", "median", "sd", "twofive_quant", "sevenfive_quant")

replicates <- 1000

continuous_variables <- c(
  'gene_gc',
  'flanking_1kb_gc',
  'lengths',
  'cm_expression',
  'ip_expression',
  'methylation',
  'H3K36me3',
  'H3K27me3',
  'H3K27ac',
  'eccdna_cov'
)

continuous_variables_subset <- c(
  'gene_gc',
  'flanking_1kb_gc',
  'lengths',
  'cm_expression',
  'ip_expression'
)

legend_labels_continuous <- c(
  'Gene GC Content',
  'Average Flanking GC Content (1kbp)',
  'Gene Length',
  'Expression: Culture (TPM)',
  'Expression: In Planta (TPM)',
  '% Methylated Cytosines',
  'Normalized H3K36me3 Signal',
  'Normalized H3K27me3 Signal',
  'Normalized H3K27ac Signal',
  'Normalized EccDNA-Seq Signal'
)

legend_labels_continuous_subset <- c(
  'Gene GC Content',
  'Average Flanking GC Content (1kbp)',
  'Gene Length',
  'Expression: Culture (TPM)',
  'Expression: In Planta (TPM)'
)

df_stats <- data.frame(matrix(ncol = length(stats)*3+2, nrow = 0))

# rbind to interleave them
colnames(df_stats) <- c(
                        "feature",
                        "host",
                        rbind(
                        paste("pav", stats_labels, sep="_"),
                        paste("conserved", stats_labels, sep="_"),
                        paste("pvalue", stats_labels, sep="_")
                        ))

for (host in c("rice", "wheat")){
  
  if (host == 'rice') {
    columns <- continuous_variables
    df_subset <- df_gene_info_rice_blast
  } else {
    columns <- continuous_variables_subset
    df_subset <- df_gene_info_wheat_blast
  }
  
  for (column in columns) {
    
    print(column)
    
    test_1 = df_subset[[column]][df_subset$lineage_pav == TRUE]
    test_2 = df_subset[[column]][df_subset$lineage_conserved == TRUE]
    
    output <- c(rbind(get_stats(test_1, stats), get_stats(test_2, stats), 
                get_all_pvalues(test_1, test_2, stats, replicates)))
    
    output <- c(column, host, output)
    
    df_stats[nrow(df_stats)+1,] <- output
  }
}

df_stats$feature <- c(legend_labels_continuous, legend_labels_continuous_subset)

write.csv(df_stats, 'per_host_continuous_param_stats.csv', row.names=FALSE)

df_gene_info_wheat_blast <-data.frame(fread('gene_info.cross_host.wheat_blast.txt', na.strings = ''))

df_gene_info_wheat_blast$host = "wheat"
  
df_gene_info_rice_blast <- data.frame(fread('gene_info.cross_host.rice_blast.txt', na.strings = ''))

df_gene_info_rice_blast$host = "rice"

# merge together after adding host names
df <- rbind(df_gene_info_rice_blast, df_gene_info_wheat_blast)

rm(df_gene_info_wheat_blast)

rm(df_gene_info_rice_blast)

df_stats <- data.frame(matrix(ncol = length(stats)*3+2, nrow = 0))

# rbind to interleave them
colnames(df_stats) <- c(
                        "feature",
                        "label",
                        rbind(
                        paste("rice", stats_labels, sep="_"),
                        paste("wheat", stats_labels, sep="_"),
                        paste("pvalue", stats_labels, sep="_")
                        ))

for (gene_status_column in c('lineage_pav', 'lineage_conserved')) {
  df_subset <- df[df[[gene_status_column]] == TRUE,]
  
  columns <- continuous_variables_subset
  
  for (column in columns) {
    
    print(column)
    
    test_1 = df_subset[[column]][df_subset$host == "rice"]
    test_2 = df_subset[[column]][df_subset$host == "wheat"]
    
    output <- c(rbind(get_stats(test_1, stats), get_stats(test_2, stats), 
                get_all_pvalues(test_1, test_2, stats, replicates)))
    
    output <- c(column, gene_status_column, output)
    
    df_stats[nrow(df_stats)+1,] <- output
  }
  
  df_stats$feature <- legend_labels_continuous_subset
  
}

write.csv(df_stats, 'per_label_continuous_param_stats.csv', row.names=FALSE)