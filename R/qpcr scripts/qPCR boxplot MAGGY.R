library(tidyverse) ## the ggplot 'universe' of packages
library(reshape2)

##import csv file, split Sample Name column, remove no template control
input_file <- "MAGGY ROS stress 2.csv" ### put the experiment name here
raw_data <- read_csv(input_file, skip=47, col_types = cols_only(
  'Sample Name' = col_character(),
  'Target Name' = col_character(),
  'CT' = col_number(),
  'Ct Mean'= col_number(),
  'Ct SD' = col_number()
))
raw_data <- as.data.frame(raw_data)
##separate here add labels as you want
raw_data <- separate(raw_data, 'Sample Name', c('treatment', 'replicate'), sep = ' ', extra = 'merge', remove = FALSE)
raw_data <- separate(raw_data, 'replicate', c('bio', 'tech'), sep = 1, extra = 'merge', remove = FALSE)

colnames(raw_data) <- c('samplename', 'treatment', 'replicate', 'bio', 'tech', 'targetname', 'ct', 'fold_change', 'ctsd')
guy11control <- subset(raw_data, treatment == 'Guy11')
raw_data <- subset(raw_data, treatment != 'Guy11')## remove controls
raw_data <- subset(raw_data, treatment != 'no')
raw_data <- subset(raw_data, select = -c(ct))
raw_data <- unique(raw_data[])
raw_data_pbr322 <- subset(raw_data, targetname == 'pBR322')
raw_data_maggy <- subset(raw_data, targetname == 'MAGGY5')

##calculate dCT
raw_data_pbr322[,7] <- raw_data_maggy[,7] - raw_data_pbr322[,7] ## dCT
raw_data_pbr322[,8] <- sqrt(raw_data_maggy[,8] ^ 2 + raw_data_pbr322[,8] ^ 2) ## get SD for dCT
raw_data_pbr322_average <- aggregate(raw_data_pbr322, by = list(raw_data_pbr322$treatment, raw_data_pbr322$bio), FUN = function (x) mean(as.numeric(as.character(x))))

control_dCT <- subset(raw_data_pbr322, treatment == 'CT')
control_dCT_average <- aggregate(control_dCT, by = list(control_dCT$treatment), FUN = function(x) mean(as.numeric(as.character(x))))
raw_data_pbr322_average[,9] <- raw_data_pbr322_average[,9] - control_dCT_average[,8]
raw_data_pbr322_average[,10] <- sqrt(raw_data_pbr322_average[,10] ^ 2 + control_dCT_average[,9] ^ 2)
raw_data_pbr322_average[,11] <- raw_data_pbr322_average[,9] + raw_data_pbr322_average[,10]
raw_data_pbr322_average[,12] <- raw_data_pbr322_average[,9] - raw_data_pbr322_average[,10]
raw_data_pbr322_average[,9] <- 2 ^ -raw_data_pbr322_average[,9]
raw_data_pbr322_average[,11] <- 2 ^ -raw_data_pbr322_average[,11]
raw_data_pbr322_average[,12] <- 2 ^ -raw_data_pbr322_average[,12]
raw_data_pbr322_average[,9] <- log(raw_data_pbr322_average[,9], 10)
raw_data_pbr322_average[,11] <- log(raw_data_pbr322_average[,11], 10)
raw_data_pbr322_average[,12] <- log(raw_data_pbr322_average[,12], 10)



pd <- position_dodge(0.5)

ggplot(raw_data_pbr322_average, aes(x = Group.1, y = fold_change, colour = Group.2)) +
  geom_point(aes(color = Group.2), position = pd, size = 5) + 
  geom_errorbar(aes(ymin = V11, ymax = V12), position = pd, width = 0.5, size = 2) +
  scale_x_discrete(limits=c('CT', 'HP', 'PQ', 'IF')) +
  theme(axis.text = element_text(size = 30)) +
  theme(legend.text=element_text(size=25), legend.title=element_text(size=25)) +
  labs(color = 'Replicate')


# position_dodge(0.2) to set jitter

# geom_boxplot(aes(fill = bio))

# size = 1.5 sets line width in geom_boxplot()

#  stat_boxplot(geom = 'errorbar', size = 1.5) +

# + geom_point(aes(color = replicate), size = 3)

ggplot(raw_data_maggy, aes(x = treatment, y = fold_change), lwd=25) +
  geom_boxplot() + geom_point(aes(color = replicate)) + labs(y = 'ct')