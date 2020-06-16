library(qqman)
setwd("C:/globus/IF_peak/")
name="RC_1AandG3.manhattanplots.table"
my_data <- read.delim(name)
tiff(paste(name, ".tiff"), units="in", width=25, height=5, res=300)
manhattan(my_data, chr = "CHROMOSOME", bp = "BASE", snp = "BASE", p = "COUNT", logp=FALSE, ylab= NA , xlab= NA , genomewideline=FALSE, suggestiveline = FALSE, main = NA, annotateTop = TRUE, cex.lab=1.5, cex.axis = 1.5, ylim = c(0,12000), col = c("blue4", "orange3"))
dev.off()
