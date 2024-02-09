#this is with thiamine
#WT and double 
directory <- '/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone/trywtanddoubleonly/'
#WT and sen1
directory <- '//Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone/trywtandsenonly/'
#WT and dbl8
directory <- '/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone/trywtanddblonly/'


#DS and dbl8!!
directory <- '/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone/tryDBLandDSonly'
#sen and dbl8!!
directory <- '/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone/trySENandDBLonly'



#this is w/o thiamine --> therefore 
#wt vs sen w/o thiamine
directory <-"/Users/patricfernandez/Documents/python/rna_seq_rep_1_and_2_no_thia/herewegoagain/allofthem/theoutput/noTHIA_wt_vs_sen/"

#wt vs dbl8 w/o thiamine
directory <-"/Users/patricfernandez/Documents/python/rna_seq_rep_1_and_2_no_thia/herewegoagain/allofthem/theoutput/noTHIA_wt_vs_dbl8/"

#wt vs DS w/o thiamine
directory <-"/Users/patricfernandez/Documents/python/rna_seq_rep_1_and_2_no_thia/herewegoagain/allofthem/theoutput/noTHIA_wt_vs_ds/"

####let's do an extra WT vs WT +thia
directory <-"/Users/patricfernandez/Documents/python/rna_seq_rep_1_and_2_no_thia/herewegoagain/allofthem/theoutput/wtnoTHIAvswtTHIA/"


sampleFiles <- grep("treated",list.files(directory),value=TRUE)
sampleCondition <- sub("(.*treated).*","\\1",sampleFiles)
#sampleReplicate <- sub(".*_(\\d+)_treated\\.txt", "\\1", sampleFiles)
sampleReplicate <- sub(".*_(\\d+)_treated\\.txt|.*_(\\d+)_untreated\\.txt", "\\1\\2", sampleFiles)
print(sampleReplicate)

newSampleNames <- c("WT", "WTa", "ds", "dsa")
newsampleReplicate <- c("wt", "wt", "ds", "ds")
#newsampleCondition <- c("one", "two", "one", "two", "one", "two", "one", "two")

sampleTable <- data.frame(sampleName = newSampleNames,
                          fileName = sampleFiles,
                          replicate = newsampleReplicate)
#sampleTable$condition <- factor(sampleTable$condition)
sampleTable$replicate <- factor(sampleTable$replicate)

print(sampleTable)

#####aoecificallyforwtvswt
sampleFiles <- grep("treated",list.files(directory),value=TRUE)
sampleCondition <- sub("(.*treated).*","\\1",sampleFiles)
#sampleReplicate <- sub(".*_(\\d+)_treated\\.txt", "\\1", sampleFiles)
sampleReplicate <- sub(".*_(\\d+)_treated\\.txt|.*_(\\d+)_untreated\\.txt", "\\1\\2", sampleFiles)
print(sampleReplicate)

newSampleNames <- c("WTthia", "WTz", "WTthiaA", "WTaZ")
newsampleReplicate <- c("WTthia", "wt", "WTthia", "wt")
#newsampleCondition <- c("one", "two", "one", "two", "one", "two", "one", "two")

sampleTable <- data.frame(sampleName = newSampleNames,
                          fileName = sampleFiles,
                          replicate = newsampleReplicate)
#sampleTable$condition <- factor(sampleTable$condition)
sampleTable$replicate <- factor(sampleTable$replicate)

print(sampleTable)
#########


library("DESeq2")
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ replicate)

dds
print(dds)


smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)
plotMA(res, ylim=c(-2,2))




resLFC <- lfcShrink(dds, coef='replicate_wt_vs_ds', type="apeglm")
results_table <- as.data.frame(resLFC)

# Filter results based on p-value threshold (e.g., 0.05)
significant_results <- subset(results_table, padj < 0.05)

library(ggplot2)

# Plot MA plot for significant results
ggplot(significant_results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  ggtitle("MA Plot of Significant Results") +
  xlab("log2 Fold Change") +
  ylab("-log10(Adjusted p-value)")

######
# Assuming 'dds' is your DESeqDataSet object

# Run DESeq analysis
dds <- DESeq(dds)

# Specify the contrast of interest
contrast <- "replicate_wt_vs_ds"

# Run the analysis and shrink log2 fold changes
resLFC <- lfcShrink(dds, coef = contrast, type = "apeglm")

# Extract results
results_table <- as.data.frame(resLFC)

# Filter results based on p-value threshold (e.g., 0.05)
significant_results <- subset(results_table, padj < 0.05)
non_significant_results <- subset(results_table, padj >= 0.05)

# Plot MA plot for both significant and non-significant results
my_plot <- ggplot() +
  geom_point(data = non_significant_results, aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.2, color = "black") +
  geom_point(data = significant_results, aes(x = log2FoldChange, y = -log10(padj)), color = "red") +
  ggtitle("MA Plot of Results") +
  xlab("log2 Fold Change") +
  ylab("-log10(Adjusted p-value)") +
  theme_minimal() +
  theme(legend.position = "none") +  # Remove legend if you don't need it
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")  # Add a horizontal line for p-value threshold

ggsave("/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone/output/wt_vs_sen_results.pdf", my_plot)

write.csv(results_table, file = "/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone/output/results_table.csv", row.names = TRUE)
write.csv(significant_results, file = "/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone/output/significant_results.csv", row.names = TRUE)
write.csv(non_significant_results, file = "/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone/output/non_significant_results.csv", row.names = TRUE)


# Assuming 'dds' is your DESeqDataSet object

# Run DESeq analysis
dds <- DESeq(dds)

# Specify the contrast of interest
contrast <- "replicate_wt_vs_ds"
contrast <- "replicate_WTthia_vs_wt"

# Run the analysis and shrink log2 fold changes
resLFC <- lfcShrink(dds, coef = contrast, type = "apeglm")

# Extract results
results_table <- as.data.frame(resLFC)

# Filter results based on p-value threshold (e.g., 0.05)
significant_results <- subset(results_table, padj < 0.05)

more_significant_results <- subset(results_table, padj < 0.00000000000000000001)
non_significant_results <- subset(results_table, padj >= 0.05)

# Plot MA plot for both significant and non-significant results

ggplot() +
  geom_point(data = non_significant_results, aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.5, color = "black") +
  geom_point(data = significant_results, aes(x = log2FoldChange, y = -log10(padj)),alpha = 0.5, color = "orange") +
  geom_point(data = more_significant_results, aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.5, color = "red") +
  geom_text(data = subset(more_significant_results, !is.na(log2FoldChange) & !is.na(padj)), 
            aes(x = log2FoldChange, y = -log10(padj), label = rownames(more_significant_results)), 
            color = "red", vjust = 1.5, hjust = 0.5, size = 3) +
  ggtitle("MA Plot of Results") +
  xlab("log2 Fold Change") +
  ylab("-log10(Adjusted p-value)") +
  theme_minimal() +
  theme(legend.position = "none") +  # Remove legend if you don't need it
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")  # Add a horizontal line for p-value threshold

write.csv(significant_results, file = "/Users/patricfernandez/Downloads/output_RNA_seq_3224/wt_vs_sen_wthia", row.names = TRUE)





######################################################################################
normalized_counts <- counts(dds)

# Subset the counts for "sen" and "sena" conditions
sen_sena_counts <- normalized_counts[, c("ds", "dsa")]

sen_sena_counts_df <- as.data.frame(sen_sena_counts)

# Add baseMean column
sen_sena_counts_df$baseMean <- rowMeans(sen_sena_counts)
write.csv(sen_sena_counts_df, file = "/Users/patricfernandez/Documents/thesis/RNA-seq\ de-seq\ outputs/mean_outputs/dbldbl.csv", row.names = TRUE)







# Step 1: Normalize counts (assuming counts are in columns 'sen' and 'sena')
sen_sena_counts_df$norm_sen <- sen_sena_counts_df$sen / sum(sen_sena_counts_df$sen) * 1e6
sen_sena_counts_df$norm_sena <- sen_sena_counts_df$sena / sum(sen_sena_counts_df$sena) * 1e6

# Step 2: Calculate TPM using the mean of normalized counts
sen_sena_counts_df$TPM <- rowMeans(sen_sena_counts_df[, c("norm_sen", "norm_sena")])


# Export the data frame to a CSV file
write.csv(sen_sena_counts_df, file = "/Users/patricfernandez/Documents/python/TPM_of_both_ds_reps.csv", row.names = TRUE)



#####trying to find the stall genes in the plot
library(dplyr)


setwd('/Users/patricfernandez/Documents/python/')
txt_df <- read.table("dbl8_stall_sites_direction.txt", header = TRUE, stringsAsFactors = FALSE)
control_df <- read.csv("new_control.csv")

##########this filters ggenes
# Assuming "more_significant_results" is your existing dataframe
# Extract the index values
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
#write.csv(filtered_txt_df, file = "missing_stall_genes_wtds.csv", row.names = TRUE)

rows <- nrow(filtresults_table)
print(rows)

ggplot() +
  geom_point(data = non_significant_results, aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.5, color = "black") +
  geom_point(data = significant_results, aes(x = log2FoldChange, y = -log10(padj)),alpha = 0.5, color = "orange") +
  geom_point(data = filtered_results, aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.9, color = "red") +
  geom_text(data = subset(filtered_results, !is.na(log2FoldChange) & !is.na(padj)), 
            aes(x = log2FoldChange, y = -log10(padj), label = rownames(filtered_results)), 
            color = "red", vjust = 1.5, hjust = 0.5, size = 3) +
  geom_point(data = filtered_results_controls, aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.9, color = "seagreen") +
  geom_text(data = subset(filtered_results_controls, !is.na(log2FoldChange) & !is.na(padj)), 
            aes(x = log2FoldChange, y = -log10(padj), label = rownames(filtered_results_controls)), 
            color = "seagreen", vjust = 1.5, hjust = 0.5, size = 3) +
  ggtitle("MA Plot of Results") +
  xlab("log2 Fold Change") +
  ylab("-log10(Adjusted p-value)") +
  theme_minimal() +
  theme(legend.position = "none") +  # Remove legend if you don't need it
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")  # Add a horizontal line for p-value threshold





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





