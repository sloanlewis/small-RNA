#The code below uses modified systemPiper for counting of reads mapped using a modified GTF.
#NB: Annotations from MiRbase are off by a 100bp - hence if the alignment has been done using a recent version of fasta, please use the ENSEMBL annotations

library(GenomicFeatures)
library(BiocParallel)
# Make sure the sqlite file is in the right destination
txdb <- loadDb("./Macaca_mulatta_smallRNA.sqlite")

eByg <- exonsBy(txdb, by="gene")
bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
multicoreParam <- MulticoreParam(workers=4); register(multicoreParam); registered()

#Counting reads - assuming strand non-specific. inter.feature=FALSE as pre- and mature miRNA positions overlap
countemiR <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=TRUE, inter.feature=FALSE, singleEnd=TRUE))
countDFmiR <- sapply(seq(along=countemiR), function(x) assays(countemiR[[x]])$counts)
rownames(countDFmiR) <- names(rowRanges(countemiR[[1]]))
colnames(countDFmiR) <- names(outpaths(args))

#Output counts and RPKM
write.table(assays(countDFmiR)$counts, "results/countDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(rpkmDFmiR, "results/rpkmDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")

#Sample outlier analysis
