library(qqman)

args = commandArgs(trailingOnly=TRUE)

my_data <- read.delim(args[1])
snps_of_interest <- readLines(args[4])

tiff(args[2], units="in", width=25, height=5, res=300)
manhattan(my_data, chr = "CHROMOSOME", bp = "BASE", snp = "SNP", p = "COUNT", highlight = snps_of_interest, highlight =  logp=FALSE, ylab= NA , xlab= NA , genomewideline=FALSE, suggestiveline = FALSE, main = NA, annotateTop = TRUE, cex.lab=1.5, cex.axis = 1.5, ylim = c(0,strtoi(args[3])), col = c("blue4", "orange3"))
dev.off()