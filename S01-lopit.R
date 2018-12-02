library("Biostrings")
library("pRolocdata")
data(hyperLOPIT2015)

## master protein identifers
up1 <- sapply(strsplit(featureNames(hyperLOPIT2015), ";"), "[[", 1)

## All uniprot
fas <- readAAStringSet("./data/uniprot_sprot.fasta.gz")
fasn <- sub("\\|.+$", "", sub("sp\\|", "", names(fas)))

## Missing identifers
up2 <- up1[!up1 %in% fasn]
writeLines(up2, "up2.txt")

## Search up2 on uniprot and download fasta
## uniprot-yourlist%3AM201812023707B37BB34871DC742D8A439176A1E10303D6F.fasta.gz
fas2 <- readAAStringSet("./data/uniprot-yourlist%3AM201812023707B37BB34871DC742D8A439176A1E10303D6F.fasta.gz")

fas <- c(fas, fas2)
names(fas) <- sub("\\|.+$", "", sub("sp\\|", "", names(fas)))

## keep features that have a fasta file
sel <- up1 %in% names(fas)

hl <- hyperLOPIT2015[sel, ]
fa <- fas[up1[sel]]

## write results
res <- data.frame(exprs(hl),
                  markers = fData(hl)$markers,
                  svm.classification = fData(hl)$svm.classification,
                  svm.score = fData(hl)$svm.score,
                  tagm.classification = fData(hl)$TAGM$tagm.mcmc.allocation,
                  tagm.probability = fData(hl)$TAGM$tagm.mcmc.probability)

write.csv(res, file = "spatprot.csv")
writeXStringSet(fa, file = "spatprot.fasta")
