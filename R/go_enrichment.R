suppressMessages(library(topGO))

args = commandArgs(trailingOnly=TRUE)

geneID2GO <-  readMappings(args[1])

geneNames <- names(geneID2GO)

myInterestingGenes <- as.character(read.csv(args[2],header=F, stringsAsFactors=FALSE, colClasses="character")$V1)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames

GOdata <- new("topGOdata", ontology = args[3], allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdata

test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)

pvalFisher <- score(resultFisher)

allRes <- GenTable(GOdata, classic = resultFisher, orderBy = "weight", ranksOf = "classic", topNodes=20)

write.table(allRes, file = args[4], sep = '\t', quote = FALSE)