library(wordcloud2)
library(data.table)
library(htmlwidgets)

args = commandArgs(trailingOnly=TRUE)

df <- data.frame(fread(args[1]))
colnames(df) <- c("word", "freq")

rownames(df) <- df[,1]

mywordcloud = wordcloud2(data = df, size = 1)
saveWidget(mywordcloud, "wordcloud.html",selfcontained = F)
webshot::webshot("wordcloud.html", args[2], vwidth = 1992, vheight = 1744, delay =10)