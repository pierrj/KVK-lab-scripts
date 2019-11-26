library(ggbiplot, lib.loc="/global/home/users/pierrj/R")
library(data.table)
mydata <- data.frame(fread("notest_transposed_final_table_for_PCA.txt"), row.names=1)
mydata <- mydata[,apply(mydata, 2, var, na.rm=TRUE) != 0]
mydata.pca <- prcomp(mydata, center = TRUE, scale. = TRUE)
mydata.biorep <- c(rep("1", 3), rep("2", 3), rep("3", 2))
tiff("PCAplot.tiff", units="in", width=5, height=5, res=300)
ggbiplot(mydata.pca, var.axes = FALSE,  groups=mydata.biorep) + geom_point(aes(colour=mydata.biorep), size=4)
dev.off()

