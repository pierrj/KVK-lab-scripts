library(tidyverse) ## the ggplot 'universe' of packages
library(reshape2)

##import csv file, split Sample Name column, remove no template control
input_file <- "Fungal DNA test infection stress 1.csv" ### put the experiment name here
raw_data <- read_csv(input_file, skip=47, col_types = cols_only(
  'Sample Name' = col_character(),
  'CT' = col_number(),
  'Ct Mean'= col_number(),
  'Ct SD' = col_number()
))
raw_data <- as.data.frame(raw_data)
##separate here add labels as you want
raw_data <- separate(raw_data, 'Sample Name', c('Type', 'Timepoint', 'Replicate'), sep = ' ', extra = 'merge', remove = FALSE)

##calculate average from control samples
fungal_controls <- subset(raw_data, Type == 'Guy11')
raw_data <- subset(raw_data, Type != 'Guy11')
fungal_average <- aggregate(fungal_controls, by = list(fungal_controls$Replicate), FUN = function(x) mean(as.numeric(as.character(x))))
fungal_average <- subset(fungal_average, select = -c(Group.1))
fungal_average[1,1] <- 'fungal_average'

##calculate dCT
polished_data <- rbind(raw_data, fungal_average)
polished_data[,5] <- polished_data[,5] - polished_data[nrow(polished_data),5] ## dCT
polished_data[,8] <- 2 ^ -polished_data[,5] ## fraction
polished_data <- subset(polished_data, Timepoint != 'NA')
colnames(polished_data) <- c('samplename', 'type', 'timepoint', 'replicate', 'ct', 'ctmean', 'ctsd', 'fraction')

ggplot(polished_data, aes(x = timepoint, y = fraction)) +
  geom_boxplot() + geom_point(aes(color = replicate))