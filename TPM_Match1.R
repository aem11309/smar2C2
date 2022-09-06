IndividualPeaks <- read.table(file = 'MaizeLeaf2hrR1_TSRPlus.ClosGeneU1K.bed', sep = '\t', header = FALSE)
GeneTotal <- read.table(file = 'MaizeLeaf2hrR1_TSSPlus.GeneSum.FINAL.bed', sep = '\t', header = FALSE)
IndividualPeaks$GeneTotal=GeneTotal$V4[match(IndividualPeaks$V12, GeneTotal$V5)]
write.csv(IndividualPeaks, file = 'MaizeLeaf2hrR1_TSRGenic.GENESUM.plus.csv')
