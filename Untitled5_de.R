#########
#wtanddouble
directory <- "/Users/patricfernandez/Documents/python/htseqoutput_29.1.24/firstattempt/rename_for_R/wt_and_double"
directory <- "/Users/patricfernandez/Documents/python/htseqoutput_29.1.24/firstattempt/rename_for_R/wt_and_dbl"
directory <- "/Users/patricfernandez/Documents/python/htseqoutput_29.1.24/firstattempt/rename_for_R/wt_and_sen"
directory<- "/Users/patricfernandez/Documents/python/htseqoutput_29.1.24/firstattempt/rename_for_R/wtonlythiaviano"

sampleFiles <- grep("treated",list.files(directory),value=TRUE)
sampleCondition <- sub("(.*treated).*","\\1",sampleFiles)
#sampleReplicate <- sub(".*_(\\d+)_treated\\.txt", "\\1", sampleFiles)
sampleReplicate <- sub(".*_(\\d+)_treated\\.txt|.*_(\\d+)_untreated\\.txt", "\\1\\2", sampleFiles)
print(sampleReplicate)

newSampleNames <- c("WT", "WTa", "ds", "dsa")
newsampleReplicate <- c("wt", "wt", "ds", "ds")
newsampleCondition <- c("one", "two", "one", "two")

sampleTable <- data.frame(sampleName = newSampleNames,
                          fileName = sampleFiles,
                          condition = newsampleCondition,
                          replicate = newsampleReplicate)
#sampleTable$condition <- factor(sampleTable$newsampleCondition)
#sampleTable$replicate <- factor(sampleTable$replicate)

print(sampleTable)


library("DESeq2")
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ condition + replicate)

dds
print(dds)


#smallestGroupSize <- 3
#keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
#dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)
plotMA(res, ylim=c(-2,2))



# Specify the contrast of interest
contrast <- "replicate_wt_vs_ds"

# Run the analysis and shrink log2 fold changes
resLFC <- lfcShrink(dds, coef = contrast, type = "apeglm")

# Extract results
results_table <- as.data.frame(resLFC)

# Filter results based on p-value threshold (e.g., 0.05)
significant_results <- subset(results_table, padj < 0.05)

more_significant_results <- subset(results_table, padj < 0.00000000000000000001)
non_significant_results <- subset(results_table, padj >= 0.05)

existing_index <- rownames(significant_results)
existing_index <- gsub("a", ".", existing_index)
# Print the modified vector
print(existing_index)
# Check for common values between "Systematic ID" column and existing index
common_values <- txt_df %>%
  filter(`ID` %in% existing_index)
# Print or do further analysis with common_values
print(common_values)



###filters dignificant_results
rownames(significant_results) <- gsub("a", ".", rownames(significant_results))
rownames(results_table) <- gsub("a", ".", rownames(results_table))

filtered_results <- significant_results[rownames(significant_results) %in% txt_df$ID , ]
filtered_results_controls <- significant_results[rownames(significant_results) %in% control_df$ID , ]


#can you try to check that all stalls are in the restults?
filtresults_table <- results_table[rownames(results_table) %in% txt_df$ID , ]
filtered_txt_df <- txt_df[!(txt_df$ID %in% rownames(results_table)), ]



ggplot() +
  geom_point(data = non_significant_results, aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.99, color = "black") +
  geom_point(data = significant_results, aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.5, color = "grey55") +
  geom_point(data = filtered_results, aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.9, color = "red") +
  geom_text(data = subset(filtered_results, !is.na(log2FoldChange) & !is.na(padj)), 
            aes(x = log2FoldChange, y = -log10(padj), label = rownames(filtered_results)), 
            color = "red", vjust = 1.5, hjust = 0.5, size = 3) +
  geom_point(data = filtered_results_controls, aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.9, color = "seagreen2") +
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
  ylim(c(-5, 80))    # Set y-axis limits


write.csv(significant_results, file = "/Users/patricfernandez/Downloads/output_RNA_seq_3224/wtthia_vs_wt_wthiaaa.csv", row.names = TRUE)






normalized_counts <- counts(dds)

# Subset the counts for "sen" and "sena" conditions
sen_sena_counts <- normalized_counts[, c("ds", "dsa")]

sen_sena_counts_df <- as.data.frame(sen_sena_counts)

# Add baseMean column
sen_sena_counts_df$baseMean <- rowMeans(sen_sena_counts)
write.csv(sen_sena_counts_df, file = "/Users/patricfernandez/Downloads/output_RNA_seq_3224/sensencounts.csv", row.names = TRUE)


