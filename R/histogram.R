args = commandArgs(trailingOnly=TRUE)

chol <- read.table(args[1], header = FALSE)
tiff(args[1], units="in", width=5, height=5, res=300)
hist(chol$V1, breaks=1000)
dev.off()