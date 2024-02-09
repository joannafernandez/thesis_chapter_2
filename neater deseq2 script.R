directory <- "/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone/pleaseplease"


directory <- "/Users/patricfernandez/Documents/python/htseqoutput_29.1.24/firstattempt/rename_for_R/with_thia_all/"
sampleFiles <- grep("treated",list.files(directory),value=TRUE)
sampleCondition <- sub("(.*treated).*","\\1",sampleFiles)
#sampleReplicate <- sub(".*_(\\d+)_treated\\.txt", "\\1", sampleFiles)
sampleReplicate <- sub(".*_(\\d+)_treated\\.txt|.*_(\\d+)_untreated\\.txt", "\\1\\2", sampleFiles)
print(sampleReplicate)

newSampleNames <- c("WT", "WTa", "sen1", "sen1a", "dbl8", "dbl8a", "ds", "dsa")
newsampleReplicate <- c("wt", "wt", "sen", "sen", "dbl", "dbl", "double", "double")
newsampleCondition <- c("one", "two", "one", "two", "one", "two", "one", "two")

sampleTable <- data.frame(sampleName = newSampleNames,
                          fileName = sampleFiles,
                          condition = newsampleCondition,
                          replicate = newsampleReplicate)
sampleTable$condition <- factor(sampleTable$condition)
sampleTable$replicate <- factor(sampleTable$replicate)

print(sampleTable)


library("DESeq2")
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ condition + replicate)


#smallestGroupSize <- 3
#keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
#dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)

####now its time to do some variance stuff

vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
#ntd <- normTransform(dds)

#meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
#meanSdPlot(assay(dds))



res

wtvssen1 <- results(dds, contrast=c("replicate","wt","sen"))
wtvsdbl <- results(dds, contrast=c("replicate","wt","dbl"))
wtvsdouble <- results(dds, contrast=c("replicate","wt","double"))


resultsNames(dds)

plotMA(wtvssen1, ylim=c(-2,2))
plotMA(wtvsdbl, ylim=c(-2,2))
plotMA(wtvsdouble, ylim=c(-2,2))

library(ggplot2)
####so this makes the shrink and volcano plot

#let's try to do this for all my contrast levels 
#this is the contrast for full data set cetween "rpelicates" --> known as condition here
resLFC <- lfcShrink(dds, coef='condition_two_vs_one', type="apeglm")
results_table <- as.data.frame(resLFC)
significant_results <- subset(results_table, padj < 0.05)



# Plot MA plot for significant results
ggplot(significant_results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  ggtitle("MA Plot of Significant Results") +
  xlab("log2 Fold Change") +
  ylab("-log10(Adjusted p-value)")



###it looks like we're gonna have to individually apply this to the folders for it to work

###here we can specify the contrast
results_shrunken <- lfcShrink(dds, contrast=c("replicate","wt","sen"))

# Display the shrunken results
print(results_shrunken)


