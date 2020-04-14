#library(ggbiplot, lib.loc="/global/home/users/pierrj/R")
library(ggbiplot)
library(data.table)
input="techrepsmerged_libnormalized"
file=paste("tableforPCA_", input, ".txt", sep="")
mydata <- data.frame(fread(file), row.names=1)
mydata <- mydata[,apply(mydata, 2, var, na.rm=TRUE) != 0]
mydata.pca <- prcomp(mydata, center = TRUE, scale. = TRUE)

##this needs to be set manually##
mydata.biorep <- c(rep("CT", 3), rep("G3", 3), rep("IF", 3), rep("PQ", 3))
##this needs to be set manually##


pc_x=1
pc_y=2
pc_choice=c(pc_x, pc_y)
tiff(paste("PCA_", "PCs", pc_x, pc_y, "_", input, ".tif", sep=""), units="in", width=5, height=5, res=300)
ggbiplot(mydata.pca, var.axes = FALSE,  groups=mydata.biorep) + geom_point(aes(colour=mydata.biorep), size=4)
dev.off()


string="G3"
#row.names.remove <- c(paste(string, "_1A", sep=""), paste(string, "_1B", sep=""), paste(string, "_1C", sep=""), paste(string, "_2A", sep=""), paste(string, "_2B", sep=""), paste(string, "_2C", sep=""), paste(string, "_3A", sep=""), paste(string, "_3B", sep=""), paste(string, "_3C", sep=""))
row.names.remove <- c(paste(string, "_1", sep=""), paste(string, "_2", sep=""), paste(string, "_3", sep=""))
mydata_removed <- mydata[!(row.names(mydata) %in% row.names.remove), ]
mydata_removed <- mydata_removed[,apply(mydata_removed, 2, var, na.rm=TRUE) != 0]
mydata_removed.pca <- prcomp(mydata_removed, center = TRUE, scale. = TRUE)

##this needs to be set manually##
mydata_removed.biorep <- c(rep("CT", 3), rep("IF", 3), rep("PQ", 3))
##this needs to be set manually##

tiff(paste("PCA_no", string ,"_", "PCs", pc_x, pc_y, "_", input, ".tif", sep=""), units="in", width=5, height=5, res=300)
ggbiplot(mydata_removed.pca, var.axes = FALSE, choices=pc_choice  ,groups=mydata_removed.biorep) + geom_point(aes(colour=mydata_removed.biorep), size=4)
dev.off()