library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

subset_file = args[1]

reference_file = args[2]

highlight_file = args[3]

plot_name = args[4]

pvaluefile = args[5]

subset <- read.table(subset_file,header=FALSE,sep='\t')
subset$V5 <- NA
subset$V6 <- NA

reference <- read.table(reference_file,header=FALSE,sep='\t')

highlight <- read.table(highlight_file, header=FALSE,sep='\t')

f_5 = function(x, output) {
  chrom = x[[1]]
  start = as.numeric(x[[2]])
  end = as.numeric(x[[3]])
  id = x[[4]]
  return(abs(max(reference$V3[reference$V1 == chrom & reference$V4 != id & reference$V3 <= start] - start)))
}

f_3 = function(x, output) {
  chrom = x[[1]]
  start = as.numeric(x[[2]])
  end = as.numeric(x[[3]])
  id = x[[4]]
  return(min(reference$V2[reference$V1 == chrom & reference$V4 != id & reference$V2 >= end] - end))
}


subset$V5 <- apply(subset, 1, f_5)

subset$V6 <- apply(subset, 1, f_3)

subset_noinf <- subset[subset$V5 != Inf & subset$V6 != Inf,]

subset_tohighlight <- subset_noinf[subset_noinf$V4 %in% highlight$V1,]

subset_noinf_log <- subset_noinf
subset_noinf_log[,5:6] <- log(subset_noinf[5:6], 10)
subset_tohighlight_log <- subset_tohighlight
subset_tohighlight_log[,5:6] <- log(subset_tohighlight[5:6], 10)

p <- ggplot(subset_noinf_log, aes(x=V5, y=V6)) + geom_point() + geom_point(data=subset_tohighlight_log, aes(x=V5, y=V6), color='red')

ggsave(plot_name, plot = p, device="pdf")

permutation_5 = replicate(10000, {
  sample_small <- subset_noinf[sample(nrow(subset_noinf), nrow(highlight) , replace = FALSE), ]
  sample_large <- subset_noinf[! rownames(subset_noinf) %in% rownames(sample_small), ]
  mean(sample_small$V5)-mean(sample_large$V5)
})

permutation_3 = replicate(10000, {
  sample_small <- subset_noinf[sample(nrow(subset_noinf), nrow(highlight) , replace = FALSE), ]
  sample_large <- subset_noinf[! rownames(subset_noinf) %in% rownames(sample_small), ]
  mean(sample_small$V6)-mean(sample_large$V6)
})

observed_5 <- mean(subset_tohighlight$V5) - mean(subset_noinf$V5)
observed_3 <- mean(subset_tohighlight$V6) - mean(subset_noinf$V6)

p.value_5 <- mean(permutation_5 > observed_5)
p.value_3 <- mean(permutation_3 > observed_3)

line_5 <- paste("the p-value for 5' distance is ", p.value_5, sep ="")
line_3 <- paste("the p-value for 3' distance is ", p.value_3, sep ="")

fileConn<-file(pvaluefile)
writeLines(c(line_5,line_3), fileConn)
close(fileConn)
