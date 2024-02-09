#directory <- "/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone/pleaseplease"


directory <- "/Users/patricfernandez/Documents/python/rna_seq_rep_1_and_2_no_thia/herewegoagain/allofthem/theoutput/"
sampleFiles <- grep("treated",list.files(directory),value=TRUE)
sampleCondition <- sub("(.*treated).*","\\1",sampleFiles)
#sampleReplicate <- sub(".*_(\\d+)_treated\\.txt", "\\1", sampleFiles)
sampleReplicate <- sub(".*_(\\d+)_treated\\.txt|.*_(\\d+)_untreated\\.txt", "\\1\\2", sampleFiles)
print(sampleReplicate)

newSampleNames <- c("WTtz", "WTz", "WTtaz", "WTaz", "sen1tz", "sen1z", "sen1taz", "sen1az", "dbl8tz", "dbl8z", "dbl8taz", "dbl8az", "dstz", "dsz", "dstaz", "dsaz")
newsampleReplicate <- c("wtt", "wt", "wtt", "wt", "sen1t", "sen1", "sen1t", "sen1", "dbl8t", "dbl8", "dbl8t", "dbl8", "dst", "ds", "dst", "ds")

newsampleReplicate <- c("wtt", "wtt", "wtt", "wtt", "sen1t", "sen1t", "sen1t", "sen1t", "dbl8t", "dbl8t", "dbl8t", "dbl8t", "dst", "dst", "dst", "dst")

newsampleCondition <- c("wthiamine", "wo", "wthiamine", "wo", "wthiamine", "wo", "wthiamine", "wo","wthiamine", "wo", "wthiamine", "wo", "wthiamine", "wo", "wthiamine", "wo")

#newsampleCondition <- c("norm", "norm", "norm", "mut", "mut", "mut", "mut", "mut","mut", "mut", "mut", "mut", "mut", "mut", "mut", "mut")

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


resultsNames(dds)
res3 <- results(dds, contrast=c("condition","wthiamine","wo"))

ntd <- normTransform(dds)

meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(dds))

vsd <- vst(dds, blind=FALSE)
plotPCA(vsd,intgroup = c("condition"))



vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

meanSdPlot(assay(vsd))

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$replicate, colData(dds)$condition, sep="-")
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "RdYlGn")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


library("pheatmap")

# Assuming colData(dds) contains a column named "replicate"
replicate_colors <- c("wtt" = "black", "sen1t" = "steelblue", "dbl8t" = "orange", 'dst' = 'tomato')
replicate_colors <- replicate_colors[as.character(colData(dds)$replicate)]

select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[, c("replicate", "condition")])

# Use the annotation_colors argument to specify the colors for the "replicate" column
pheatmap(
  assay(dds)[select,],
  cluster_rows = FALSE,
  show_rownames = FALSE,
  cluster_cols = FALSE,
  annotation_col = df,
  annotation_colors = list(replicate = replicate_colors),
  color = colorRampPalette(brewer.pal(9, "RdYlBl"))(100)  # Change the color scale as needed
)






plotPCA(vsd,intgroup = c("replicate"))

ntdresLFC <- lfcShrink(dds, coef='condition_wthiamine_vs_wo', type="apeglm")
results_table <- as.data.frame(ntdresLFC)

# Filter results based on p-value threshold (e.g., 0.05)
significant_results <- subset(results_table, padj < 0.05)

library(ggplot2)

# Plot MA plot for significant results
ggplot(significant_results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  ggtitle("MA Plot of Significant Results") +
  xlab("log2 Fold Change") +
  ylab("-log10(Adjusted p-value)")


# Specify the contrast of interest
contrast <- "condition_wthiamine_vs_wo"

# Run the analysis and shrink log2 fold changes
resLFC <- lfcShrink(dds, coef = contrast, type = "apeglm")

# Extract results
results_table <- as.data.frame(resLFC)

# Filter results based on p-value threshold (e.g., 0.05)
significant_results <- subset(results_table, padj < 0.05)

more_significant_results <- subset(results_table, padj < 0.00000000000000000001)
non_significant_results <- subset(results_table, padj >= 0.05)

ggplot() +
  geom_point(data = non_significant_results, aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.99, color = "black") +
  geom_point(data = significant_results, aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.5, color = "grey55") +

  #geom_text(data = subset(filtered_results_controls, !is.na(log2FoldChange) & !is.na(padj)), 
  #         aes(x = log2FoldChange, y = -log10(padj), label = rownames(filtered_results_controls)), 
  #        color = "seagreen", vjust = 1.5, hjust = 0.5, size = 3) +
  ggtitle("MA Plot of Results") +
  xlab("log2 Fold Change") +
  ylab("-log10(Adjusted p-value)") +
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend if you don't need it
        axis.text = element_text(size = 17),  # Set font size for axis labels
  ) +  # Remove gridlines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # Add a horizontal line for p-value threshold
  xlim(c(-6, 6)) +  # Set x-axis limits
  ylim(c(-5, 250))    # Set y-axis limits


write.csv(significant_results, file = "/Users/patricfernandez/Downloads/output_RNA_seq_3224/wthia_vs_wo.csv", row.names = TRUE)


