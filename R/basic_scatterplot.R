## USAGE ##
## options in order ##
# input file name
# x column number
# y column number
# xlab
# ylab
# plot title
# file title
# file x dim
# file y dim

library(data.table)

args = commandArgs(trailingOnly=TRUE)

df <- data.frame(fread(args[1]))

plot_title = paste(args[7], ".pdf", sep = "")

pdf(file = plot_title, width = as.numeric(args[8]), height = as.numeric(args[9]))

plot(x = df[,as.numeric(args[2])], y = df[,as.numeric(args[3])],
     xlab = args[4], ylab = args[5],
     main = args[6])

dev.off()