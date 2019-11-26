library(tidyverse) ## the ggplot 'universe' of packages
library(reshape2)

##import csv file, split Sample Name column, remove no template control
input_file <- "experiment 3 QC r script practice/data.csv" ### put the experiment name here
raw_data <- read_csv(input_file, skip=47, col_types = cols_only(
  'Sample Name' = col_character(),
  'CT' = col_number(),
  'Ct Mean'= col_number(),
  'Ct SD' = col_number()
))
raw_data <- as.data.frame(raw_data)
raw_data <- separate(raw_data, 'Sample Name', c('Genotype', 'Method'), sep = ' ', extra = 'merge', remove = FALSE)
raw_data <- subset(raw_data, Genotype != 'no')

##calculate average from gDNA samples
gDNA_data <- subset(raw_data, Method == 'gDNA')
raw_data <- subset(raw_data, Method != 'gDNA')
gDNA_average <- aggregate(gDNA_data, by = list(gDNA_data$Method), FUN = function(x) mean(as.numeric(as.character(x))))
gDNA_average <- subset(gDNA_average, select = -c(Group.1))
gDNA_average[1,1] <- 'gDNA_average'

##calculate dCT and remaining linear fraction
polished_data <- rbind(raw_data, gDNA_average)
polished_data[,4] <- polished_data[,4] - polished_data[nrow(polished_data),4] ## dCT
polished_data[,7] <- 2 ^ -polished_data[,4] ## fraction
polished_data <- subset(polished_data, Method != 'NA')
colnames(polished_data) <- c('samplename', 'genotype', 'method', 'ct', 'ctmean', 'ctsd', 'fraction')

ggplot(polished_data, aes(x = method, y = fraction)) +
  geom_boxplot() + geom_point(aes(color = samplename))

