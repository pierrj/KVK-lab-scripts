my_data <- read.delim("all_mapcounts.txt")
counts <- table(my_data$COUNT)
barplot(counts, main="mapping distribution", 
        xlab="Number of times reads map to genome")
