#library(ggbiplot, lib.loc="/global/home/users/pierrj/R")
library(ggbiplot)
library(data.table)
setwd("C:/globus/merged_tests/")
input="techrepsmerged_spikenormalized"
file=paste("tableforPCA_", input, ".txt", sep="")
mydata <- data.frame(fread(file), row.names=1)
mydata <- mydata[,apply(mydata, 2, var, na.rm=TRUE) != 0]
mydata.pca <- prcomp(mydata, center = TRUE, scale. = TRUE)

##this needs to be set manually##
 mydata.biorep <- c(rep("CT_1", 3), rep("CT_2", 3), rep("CT_3", 3), , rep("IF_1", 3), rep("IF_2", 3), rep("IF_3", 3), rep("PQ_1", 3), rep("PQ_2", 3), rep("PQ_3", 3))
# mydata.biorep <- c(rep("CT", 3), rep("IF", 3), rep("PQ", 3))

##this needs to be set manually##


pc_x=1
pc_y=2
pc_choice=c(pc_x, pc_y)
tiff(paste("PCA_", "PCs", pc_x, pc_y, "_", input, ".tif", sep=""), units="in", width=5, height=5, res=300)
ggbiplot(mydata.pca, var.axes = FALSE, choices=pc_choice,  groups=mydata.biorep) + geom_point(aes(colour=mydata.biorep), size=4)
dev.off()


string=""
string2=""
string3=""
string4=""
string5=""
string6=""
# row.names.remove <- c(paste(string, "_1A", sep=""), paste(string, "_1B", sep=""), paste(string, "_1C", sep=""), paste(string, "_2A", sep=""), paste(string, "_2B", sep=""), paste(string, "_2C", sep=""), paste(string, "_3A", sep=""), paste(string, "_3B", sep=""), paste(string, "_3C", sep=""),
#                      paste(string2, "_1A", sep=""), paste(string2, "_1B", sep=""), paste(string2, "_1C", sep=""), paste(string2, "_2A", sep=""), paste(string2, "_2B", sep=""), paste(string2, "_2C", sep=""), paste(string2, "_3A", sep=""), paste(string2, "_3B", sep=""), paste(string2, "_3C", sep=""),
#                      paste(string3, "_1A", sep=""), paste(string3, "_1B", sep=""), paste(string3, "_1C", sep=""), paste(string3, "_2A", sep=""), paste(string3, "_2B", sep=""), paste(string3, "_2C", sep=""), paste(string3, "_3A", sep=""), paste(string3, "_3B", sep=""), paste(string3, "_3C", sep=""),
#                      paste(string4, "_1A", sep=""), paste(string4, "_1B", sep=""), paste(string4, "_1C", sep=""), paste(string4, "_2A", sep=""), paste(string4, "_2B", sep=""), paste(string4, "_2C", sep=""), paste(string4, "_3A", sep=""), paste(string4, "_3B", sep=""), paste(string4, "_3C", sep=""),
#                      string5, string6)

row.names.remove <- c(paste(string, "_1", sep=""), paste(string, "_2", sep=""), paste(string, "_3", sep=""),
                      paste(string2, "_1", sep=""), paste(string2, "_2", sep=""), paste(string2, "_3", sep=""),
                      string5)
#row.names.remove <- c(string)
mydata_removed <- mydata[!(row.names(mydata) %in% row.names.remove), ]
mydata_removed <- mydata_removed[,apply(mydata_removed, 2, var, na.rm=TRUE) != 0]
mydata_removed.pca <- prcomp(mydata_removed, center = TRUE, scale. = TRUE)
mydata_removed.pca$rotation[,"PC1"]

##this needs to be set manually##
mydata_removed.biorep <- c(rep("Control", 3),rep("Infected", 3), rep("ROS", 3))
##this needs to be set manually##

tiff(paste("PCA_no", string, string2 , string3, string4, string5, string6, "_", "PCs", pc_x, pc_y, "_", input, ".tif", sep=""), units="in", width=5, height=5, res=300)
ggbiplot(mydata_removed.pca, var.axes = FALSE, choices=pc_choice  ,groups=mydata_removed.biorep ) + geom_point(aes(colour=mydata_removed.biorep), size=4) + labs(color="Treatment")
dev.off()


matrix<-as.matrix(vegdist(mydata_removed, method="euclidian"))
df <- data.frame(fread("treatments.txt"), row.names=1)
adonis(matrix ~ treatment, data = df)

