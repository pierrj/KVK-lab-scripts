library(tidyverse) ## the ggplot 'universe' of packages
library(reshape2)

##import csv file, split Sample Name column, remove no template control
input_file <- "digestion time course 1 r script practice/data.csv" ### put the experiment name here
raw_data <- read_csv(input_file, skip=47, col_types = cols_only(
  'Sample Name' = col_character(),
  'CT' = col_number(),
  'Ct Mean'= col_number(),
  'Ct SD' = col_number(),
  'Target Name' = col_character()
))
raw_data <- as.data.frame(raw_data)
raw_data <- separate(raw_data, 'Sample Name', c('Genotype', 'Method'), sep = ' ', extra = 'merge', remove = FALSE)
raw_data <- separate(raw_data, 'Target Name', c('Target', 'Name'), sep = ' ', extra = 'merge', remove = FALSE)
raw_data <- subset(raw_data, Genotype != 'no')
raw_data <- subset(raw_data, Target != 'YFP')

##calculate average from "pure" samples
pure_data <- subset(raw_data, Method == 'pure')
raw_data <- subset(raw_data, Method != 'pure')
raw_data <- raw_data[grep("dig", raw_data$Method),]
pure_average <- aggregate(pure_data, by = list(pure_data$Genotype), FUN = function(x) mean(as.numeric(as.character(x))))
pure_average[,3] <- pure_average[,1]
pure_average <- subset(pure_average, select = -c(Group.1))


##calculate dCT and remaining linear fraction
raw_data$CT <- raw_data$CT - pure_average$CT[match(raw_data$Genotype, pure_average$Genotype)]
raw_data[,10] <- 2 ^ -raw_data[,7] ## fraction

##set better column names, then assign time points to each genotype/method combination
raw_data[,11] <- NA
colnames(raw_data) <- c('samplename', 'genotype', 'method', 'target name', 'target', 'name', 'dct', 'ctmean', 'ctsd', 'fraction', 'timepoint')
for (row in 1:nrow(raw_data)) {
  genotype <- raw_data[row, 'genotype']
  method <- raw_data[row, 'method']
  if (genotype=='KA' | genotype=='GA')
  {
  if (method == 'dig 1') {
    raw_data[row,11] <- '48h'}
  else if(method == 'dig 2') {
    raw_data[row,11] <- '72h'}
  else if(method == 'dig 3') {
    raw_data[row,11] <- '96h'}
  else {next}
  }
  else
  {
    if (method == 'dig 1') {
      raw_data[row,11] <- '24h'}
    else if(method == 'dig 2') {
      raw_data[row,11] <- '48h'}
    else if(method == 'dig 3') {
      raw_data[row,11] <- '72h'}
    else {next}
  }
  }

##plot
ggplot(raw_data, aes(x = timepoint, y = fraction, colour = timepoint)) +
  geom_boxplot() + geom_point()
