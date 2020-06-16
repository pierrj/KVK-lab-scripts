library(wordcloud2)
library(data.table)
library(htmlwidgets)

Pfamtop50 <- data.frame(fread("C:/globus/pfam/G3.hmmscan.domainlist.trimmed.counted.txt"))
colnames(Pfamtop50) <- c("freq", "word")

Pfamtop50 <- Pfamtop50[c("word", "freq")]
rownames(Pfamtop50) <- Pfamtop50[,1]


mywordcloud = wordcloud2(data = Pfamtop50, size = 1)
saveWidget(mywordcloud, "test1.html",selfcontained = F)
webshot::webshot("test1.html","test1.png",vwidth = 1992, vheight = 1744, delay =10)
