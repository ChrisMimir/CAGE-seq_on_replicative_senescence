#Extended Data Figure 1a
# Step 1: Extract TPM and organize
tpm_matrix <- assay(RSE, "TPM")

anno_df <- as.data.frame(rowData(RSE)) %>%
  rownames_to_column("cluster_id") %>%
  select(cluster_id, txType)

sample_info <- as.data.frame(colData(RSE)) %>%
  rownames_to_column("sample_name") %>%
  select(sample_name, Class)

tpm_long <- as.data.frame(tpm_matrix) %>%
  rownames_to_column("cluster_id") %>%
  pivot_longer(-cluster_id, names_to = "sample_name", values_to = "TPM") %>%
  left_join(anno_df, by = "cluster_id") %>%
  left_join(sample_info, by = "sample_name") %>%
  mutate(logTPM = log2(TPM + 1)) %>%
  filter(!is.na(txType))

# Step 2: Calculate mean expression (by txType and Class)
tpm_summary <- tpm_long %>%
  group_by(txType, Class) %>%
  summarise(mean_logTPM = mean(logTPM, na.rm = TRUE), .groups = "drop")

# Optional: Ensure Class order is correct
tpm_summary$Class <- factor(tpm_summary$Class,
                            levels = c("p03", "p09", "p12", "p14", "p16", "p18"))

ggplot(tpm_summary, aes(x = Class, y = mean_logTPM, fill = Class)) +
  geom_col(color = "black", width = 0.7, alpha = 0.9) +
  facet_wrap(~ txType, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(
    "p03" = "#F8766D",
    "p09" = "#C49A00",
    "p12" = "#53B400",
    "p14" = "#00C1AA",
    "p16" = "#00B0F6",
    "p18" = "#A58AFF"
  )) +
  labs(x= "passage", y = expression("log2(TPM)")) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 13),
    panel.spacing = unit(1, "lines"),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )

ggplot(tpm_summary, aes(x = Class, y = mean_logTPM, fill = Class)) +
  geom_col(width = 0.7, alpha = 0.9, 
           # 1. Fix: Remove black border
           color = NA) +
  facet_wrap(~ txType, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(
    "p03" = "#F8766D",
    "p09" = "#C49A00",
    "p12" = "#53B400",
    "p14" = "#00C1AA",
    "p16" = "#00B0F6",
    "p18" = "#A58AFF"
  )) +
  labs(x= "passage", y = expression("log2(TPM)")) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(face = "bold"),
    # 2. Fix: Remove facet title background/border
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    panel.spacing = unit(1, "lines"),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )



#Extended Data Figure 1b
#Prepare TSS-TE data
TE <- read.table("C:/Users/Jiaqi/Documents/Code/R/aging/human_TE_TSS.txt",
                 header=TRUE, stringsAsFactors=FALSE, fill=TRUE, sep='\t')
head(TE)
TE$strand <- ifelse(TE$strand == "", "*", TE$strand)
TE_DE <- GRanges(
  seqnames = TE$chr,
  ranges = IRanges(start = TE$start, end = TE$end),
  strand = TE$strand,
  type = TE$type,
  transposon  = TE$transposon,
  gene_id = TE$gene_id,
  gene_name = TE$gene_name)
head(TE_DE)

TERSE <- RSE

normalized_counts <- counts(dds, normalized = TRUE)


assay(TERSE, "size_factor_normalized_counts") <- normalized_counts

# Extract normalized count matrix
expr_matrix <- assay(TERSE, "size_factor_normalized_counts")  # Replace with appropriate assay name

# Extract group information
group_info <- colData(TERSE)$Class  # Group information
group_info <- factor(group_info)
average_exp <- apply(expr_matrix, 1, function(x) tapply(x, group_info, mean, na.rm = TRUE))

# Convert average expression matrix to data frame
average_exp_df <- as.data.frame(t(average_exp)) 
# Add group names to column names
colnames(average_exp_df) <- levels(group_info)
# Add average expression to rowData
rowData(TERSE)$average_exp <- average_exp_df

rowData(TERSE)$type <- "none"
rowData(TERSE)$transposon <- "none"
overlaps <- findOverlaps(rowRanges(RSE), TE_DE)

# Extract overlapping row indices
TERSE_hits <- queryHits(overlaps)
TE_DE_hits <- subjectHits(overlaps)

# Merge matching annotations (handle duplicate matches)
aggregated_type <- tapply(TE_DE$type[TE_DE_hits], TERSE_hits, paste, collapse = ",")
aggregated_transposon <- tapply(TE_DE$transposon[TE_DE_hits], TERSE_hits, paste, collapse = ",")

# Align annotation information to TERSE rows
aligned_type <- rep("none", nrow(TERSE))
aligned_transposon <- rep("none", nrow(TERSE))

# Update annotation information by index
aligned_type[as.numeric(names(aggregated_type))] <- aggregated_type
aligned_transposon[as.numeric(names(aggregated_transposon))] <- aggregated_transposon

# Add aligned annotations to rowData
rowData(TERSE)$type <- aligned_type
rowData(TERSE)$transposon <- aligned_transposon

# Check results
head(rowData(TERSE))

TERSE <- TERSE[rowData(TERSE)$type != "none", ]

# Extract normalized count matrix
expr_matrix <- assay(TERSE, "size_factor_normalized_counts")  # Replace with your normalized count matrix

# Extract group information
group_info <- colData(TERSE)$Class  # Replace with your group column name

# Convert group information to factor
group_info <- factor(group_info)

# Calculate total expression sum per group
total_expression <- colSums(expr_matrix)

# Calculate total counts by group
group_totals <- tapply(total_expression, group_info, sum)

# Convert to data frame format
group_totals_df <- data.frame(
  Class = names(group_totals),
  TotalExpression = as.numeric(group_totals)
)

ggplot(group_totals_df, aes(x = Class, y = TotalExpression)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal(base_size = 15) +
  labs(
    title = "Overall Transposable elements expression",
    x = "Class",
    y = "Total Expression"
  )

TEA <- TERSE
rowData(TEA)$transposon
te_blind <- DESeqDataSet(TEA, design = ~ 1)
# Use varianceStabilizingTransformation
te_blind <- varianceStabilizingTransformation(te_blind, blind = TRUE)

plotPCA(te_blind, "Class")

TEA_seqs <- TEA %>% 
  rowRanges() %>%
  swapRanges() %>%
  unstrand() %>%
  add(500) %>% 
  getSeq(bsg, names = .)
transposon_data <- rowData(TEA)$transposon

# Split with strsplit, and fill missing values with NA
transposon_split <- lapply(strsplit(transposon_data, "\\|"), function(x) {
  if (length(x) < 3) {
    c(x, rep(NA, 3 - length(x)))  # If length < 3, pad with NA
  } else {
    x[1:3]  # If length >= 3, take only first 3 columns
  }
})

transposon_split <- do.call(rbind, transposon_split)
# Name the three columns
colnames(transposon_split) <- c("class", "family", "subfamily")

# Store the matrix directly in the transposon field
rowData(TEA)$transposon <- transposon_split
rowData(TEA)$transposon

# Get all transposon classes
transposon_classes <- unique(rowData(TEA)$transposon[, "class"])

calculate_transposon_matrix <- function(expr_matrix, group_info, transposon_classes, transposon_annotation) {
  # Initialize result matrix
  result_matrix <- matrix(0, nrow = length(transposon_classes), ncol = length(levels(group_info)))
  rownames(result_matrix) <- transposon_classes
  colnames(result_matrix) <- levels(group_info)
  
  # Calculate weighted expression for each class by group
  for (group in levels(group_info)) {
    # Filter samples for current group
    group_samples <- which(group_info == group)
    
    # Expression matrix for current group
    group_expr <- expr_matrix[, group_samples]
    
    # Calculate total expression for each transposon class
    for (class in transposon_classes) {
      # Filter TSS for current class
      class_tss <- which(transposon_annotation == class)
      
      # Calculate total weighted expression sum
      result_matrix[class, group] <- sum(group_expr[class_tss, ], na.rm = TRUE)
    }
  }
  
  return(result_matrix)
}

# Extract necessary data
expr_matrix <- assay(TEA, "size_factor_normalized_counts")  # Replace with normalized count matrix
group_info <- factor(colData(TEA)$Class)  # Replace with group information
transposon_annotation <- rowData(TEA)$transposon[, "class"]  # Transposon class annotation

# Calculate result matrix
transposon_matrix <- calculate_transposon_matrix(expr_matrix, group_info, transposon_classes, transposon_annotation)

head(transposon_matrix)

transposon_df <- as.data.frame(transposon_matrix)
transposon_df <- transposon_df %>%
  rownames_to_column(var = "class") %>%
  pivot_longer(cols = -class, names_to = "Group", values_to = "Expression")

# Create stacked bar plot
ggplot(transposon_df, aes(x = Group, y = Expression, fill = class)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 15) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Transposon Class Expression",
    x = "Passage",
    y = "Expression (TPM)",
    fill = "Transposon Class"
  )

#-------------------------------
# Extract TPM matrix (or normalized counts)
tpm_matrix <- assay(TERSE, "TPM")  # or size_factor_normalized_counts

# Extract sample information
sample_info <- as.data.frame(colData(TERSE)) %>%
  rownames_to_column("sample_name") %>%
  select(sample_name, Class)

# Extract annotation information (transposon and genomic location type)
anno_df <- as.data.frame(rowData(TERSE)) %>%
  rownames_to_column("tss_id") %>%
  mutate(
    TE_class = sapply(strsplit(transposon, "\\|"), `[`, 1),
    TE_class = ifelse(is.na(TE_class) | TE_class == "", "none", TE_class)
  )

tpm_long <- as.data.frame(tpm_matrix) %>%
  rownames_to_column("tss_id") %>%
  pivot_longer(-tss_id, names_to = "sample_name", values_to = "TPM") %>%
  left_join(anno_df, by = "tss_id") %>%
  left_join(sample_info, by = "sample_name") %>%
  filter(TE_class != "none" & !is.na(txType) & !is.na(Class)) %>%
  mutate(logTPM = log2(TPM + 1))

te_summary <- tpm_long %>%
  group_by(txType, Class) %>%
  summarise(mean_logTPM = mean(logTPM, na.rm = TRUE),
            .groups = "drop")

# Ensure Class order is consistent
te_summary$Class <- factor(te_summary$Class, 
                           levels = c("p03", "p09", "p12", "p14", "p16", "p18"))

# PCA plot consistent color scheme
my_cols <- c("p03" = "#F8766D","p09" = "#C49A00","p12" = "#53B400",
             "p14" = "#00C1AA","p16" = "#00B0F6","p18" = "#A58AFF")

ggplot(te_summary, aes(x = Class, y = mean_logTPM, fill = Class)) +
  geom_col(color = "black", width = 0.7, alpha = 0.9) +
  facet_wrap(~ txType, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = my_cols) +
  labs(
    title = "Transposable elements derived TSS by genomic location",
    x = "Cell passage",
    y = expression("mean log2(TPM)")
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )

#----------Extended Data Figure 4b--------------------------------------#
# Assuming RSE contains TSS-level reads
RSE_tss <- subset(RSE, clusterType == "TSS")

# Design matrix: time groups (e.g., 6 passages)
colData(RSE_tss)$group <- factor(colData(RSE_tss)$Class,
                                 levels = c("p03","p09","p12","p14","p16","p18"))

dds <- DESeqDataSet(RSE_tss, design = ~ group)

# Filter low-expression TSS
dds <- dds[rowSums(counts(dds)) > 10, ]

# DESeq normalization and differential analysis
dds <- DESeq(dds)

# Example: p03 vs p18
res <- results(dds, contrast = c("group", "p18", "p03"))
res <- lfcShrink(dds, contrast = c("group", "p18", "p03"), res = res, type = "ashr")

# Add symbol annotation
res$symbol <- rowData(RSE_tss)$symbol[match(rownames(res), rownames(RSE_tss))]

# Add significance markers
res$DE <- ifelse(res$padj < 0.05 & res$log2FoldChange > 1, "Up",
                 ifelse(res$padj < 0.05 & res$log2FoldChange < -1, "Down", "NS"))

# Extract all comparison results
comparisons <- combn(levels(colData(RSE_tss)$group), 2, simplify = FALSE)

# Summarize differential results for each group
summary_list <- lapply(comparisons, function(comp) {
  res <- results(dds, contrast = c("group", comp[2], comp[1]))
  res$comparison <- paste(comp[2], "vs", comp[1])
  res$DE <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1,
                   ifelse(res$log2FoldChange > 0, "Up", "Down"), "NS")
  res
})



summary_list <- lapply(comparisons, function(comp) {
  res <- results(dds, contrast = c("group", comp[2], comp[1]))
  res <- lfcShrink(dds, contrast = c("group", comp[2], comp[1]), res = res, type = "ashr")
  
  res_df <- as.data.frame(res)   # ✅ Convert to data.frame
  res_df$comparison <- paste(comp[2], "vs", comp[1])
  res_df$gene_id <- rownames(res_df)
  res_df$DE <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                      ifelse(res_df$log2FoldChange > 0, "Up", "Down"), "NS")
  return(res_df)
})

dt <- do.call(rbind, summary_list)

# Global statistics
summary_table <- dt %>%
  group_by(comparison, DE) %>%
  summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = DE, values_from = n, values_fill = 0)

knitr::kable(summary_table,
             caption = "Global summary of differentially expressed TSSs.")

links <- cor_links

sig_tss_by_group <- dt %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  group_by(comparison) %>%
  summarise(sig_tss = list(unique(gene_id)))

tss_ids <- names(RSE_tss)

link_pairs <- purrr::map2_df(
  names(links),
  links$revmap,
  ~ data.frame(
    enhancer = .x,
    TSS = tss_ids[.y]
  )
)

enh_summary <- sig_tss_by_group %>%
  mutate(
    enh_count = purrr::map_int(sig_tss, ~ {
      sum(link_pairs$TSS %in% .x)
    })
  ) %>%
  select(comparison, enh_count)

ggplot(enh_summary, aes(x = comparison, y = enh_count)) +
  geom_col(fill = "#009E73", width = 0.7) +
  theme_bw(base_size = 14) +
  labs(
    x = "Comparison",
    y = "Number of enhancers linked to DE TSSs"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

selected_comparisons <- c("p09 vs p03", "p12 vs p09", "p14 vs p12", "p16 vs p14", "p18 vs p16")

# Extract and sort
enh_plot_df <- enh_summary %>%
  filter(comparison %in% selected_comparisons) %>%
  mutate(comparison = factor(comparison, levels = selected_comparisons))

# Custom colors
my_colors <- c(
  "#F8766D",
  "#C49A00",
  "#53B400",
  "#00C1AA",
  "#00B0F6"
)

# Create plot
p <- ggplot(enh_plot_df, aes(x = comparison, y = enh_count, fill = comparison)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = my_colors) +
  labs(
    x = "Contrast",
    y = "Number of enhancers linked to DE TSSs"
  ) +
  theme_bw(base_size = 14) +
  theme(
    # Remove right and top borders
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    # Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Tilt X-axis text
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Remove legend
    legend.position = "none"
  )

# Display
print(p)
