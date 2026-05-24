library(dplyr)
library(scales)
library(cmapR)
library(GEOquery)
library(edgeR)
library(broom)
library(limma)
library(ggplot2)
# Keep only for publication-level plotting

# Assuming data is already loaded: tissue_exp, sample_annot, subject_annot

expr_mat <- tissue_exp@mat    
row_ann  <- tissue_exp@rid      
col_ann  <- tissue_exp@cid

# Target gene
target_gene_id <- "ENSG00000123095.6"
Target_gene <- "BHLHE41"

#=======================================================
# 2. Build metadata (Base R syntax)
#=======================================================

# Extract and rename sample_info
# Note: Strongly recommend adding ischemic time SMTSISCH as a covariate
sample_info <- sample_annot[, c("SAMPID", "SMTSD", "SMRIN", "SMTSISCH")]
colnames(sample_info) <- c("SAMPID", "Tissue", "RIN", "IschemicTime")
sample_info$SUBJID <- substr(sample_info$SAMPID, 1, 10)

# Extract and rename subject_info
subject_info <- subject_annot[, c("SUBJID", "AGE", "SEX")]
colnames(subject_info) <- c("SUBJID", "Age", "Sex")

# Merge metadata and filter samples present in column annotations
meta_full <- merge(sample_info, subject_info, by = "SUBJID")
meta_full <- meta_full[meta_full$SAMPID %in% col_ann, ]

# Age mapping
age_map <- c(
  "20-29" = 25, "30-39" = 35, "40-49" = 45,
  "50-59" = 55, "60-69" = 65, "70-79" = 75
)
meta_full$Age_num <- age_map[as.character(meta_full$Age)]

# Clean expression matrix
expr_mat <- as.matrix(expr_mat)
expr_mat <- expr_mat[!is.na(rownames(expr_mat)), !is.na(colnames(expr_mat))]

#=======================================================
# 3. Tissue-by-tissue analysis (with QC and covariate adjustment)
#=======================================================

results_list <- list()
plot_data_list <- list()

tissues <- unique(meta_full$Tissue)

for (tissue in tissues) {
  cat("Processing:", tissue, "\n")
  
  #----- QC -----
  is_brain <- grepl("^Brain-", tissue)
  
  # Extract samples for this tissue, excluding missing age or RIN values
  meta_t <- meta_full[meta_full$Tissue == tissue & 
                        !is.na(meta_full$Age_num) & 
                        !is.na(meta_full$RIN), ]
  
  if (is_brain) {
    meta_t <- meta_t[meta_t$RIN >= 3, ]
  } else {
    meta_t <- meta_t[meta_t$RIN >= 6, ]
  }
  
  if (nrow(meta_t) < 60) {
    cat("  -> Skip (QC sample < 60):", nrow(meta_t), "\n")
    next
  }
  
  #----- Extract expression matrix and filter -----
  expr_t <- expr_mat[, meta_t$SAMPID, drop = FALSE]
  keep_genes <- rowSums(expr_t >= 1) >= 5
  expr_t_filtered <- expr_t[keep_genes, , drop = FALSE]
  
  if (nrow(expr_t_filtered) < 100) { next }
  
  #----- TMM & Voom -----
  dge <- DGEList(counts = expr_t_filtered)
  dge <- calcNormFactors(dge, method = "TMM")
  
  design_voom <- model.matrix(~ 1, data = meta_t)
  v <- voom(dge, design_voom, plot = FALSE)
  
  if (!target_gene_id %in% rownames(v$E)) { next }
  
  #----- Extract target gene expression data -----
  target_voom <- v$E[target_gene_id, ]
  
  df_plot <- meta_t
  # Use match to ensure sample alignment
  df_plot$target_norm <- target_voom[match(df_plot$SAMPID, colnames(v$E))]
  
  #----- Multiple linear regression (Top Journal Standard) -----
  # Dynamically determine if sex can be used as a covariate (exclude ovary, prostate, etc.)
  has_sex_covariate <- length(unique(df_plot$Sex)) > 1
  
  if (has_sex_covariate) {
    fit <- lm(target_norm ~ Age_num + Sex + RIN, data = df_plot) 
  } else {
    fit <- lm(target_norm ~ Age_num + RIN, data = df_plot)
  }
  
  # Extract Age_num statistics from model results
  coef_summary <- summary(fit)$coefficients
  if ("Age_num" %in% rownames(coef_summary)) {
    res_df <- data.frame(
      Tissue = tissue,
      estimate = coef_summary["Age_num", "Estimate"],
      p.value = coef_summary["Age_num", "Pr(>|t|)"],
      N = nrow(df_plot),
      stringsAsFactors = FALSE
    )
    results_list[[tissue]] <- res_df
  }
  
  plot_data_list[[tissue]] <- df_plot
}

#=======================================================
# 4. Merge results and FDR correction (Base R implementation)
#=======================================================

if (length(results_list) == 0) stop("No tissues passed QC.")

# Merge list into data frame
lm_results <- do.call(rbind, results_list)
plot_df <- do.call(rbind, plot_data_list)

# Calculate FDR and add significance stars
lm_results$FDR <- p.adjust(lm_results$p.value, method = "BH")
lm_results$sig_star <- ifelse(lm_results$FDR < 0.001, "***",
                              ifelse(lm_results$FDR < 0.01, "**",
                                     ifelse(lm_results$FDR < 0.05, "*", "")))

# Build plot text labels
lm_results$label <- paste0("N = ", lm_results$N, 
                           "\nβ = ", signif(lm_results$estimate, 3), 
                           "\nFDR = ", signif(lm_results$FDR, 3), " ", lm_results$sig_star)

# Dynamically calculate x, y coordinates for labels to prevent occlusion (using Base R aggregate)
y_pos <- aggregate(target_norm ~ Tissue, data = plot_df, FUN = function(x) max(x, na.rm = TRUE) * 1.05)
colnames(y_pos) <- "y_text"

x_pos <- aggregate(Age_num ~ Tissue, data = plot_df, FUN = function(x) min(x, na.rm = TRUE))
colnames(x_pos) <- "x_text"

# Merge coordinates into lm_results
# 1. Calculate y coordinates and force column names
y_pos <- aggregate(target_norm ~ Tissue, data = plot_df, FUN = function(x) max(x, na.rm = TRUE) * 1.05)
colnames(y_pos) <- c("Tissue", "y_text") # Force first column to be Tissue, second to be y_text

# 2. Calculate x coordinates and force column names
x_pos <- aggregate(Age_num ~ Tissue, data = plot_df, FUN = function(x) min(x, na.rm = TRUE))
colnames(x_pos) <- c("Tissue", "x_text") # Force first column to be Tissue, second to be x_text

# 3. Merge x and y coordinates
pos_df <- merge(x_pos, y_pos, by = "Tissue")

# 4. Merge coordinates back to statistical results (check if lm_results also has Tissue column, which it normally does)
lm_results <- merge(lm_results, pos_df, by = "Tissue")

#=======================================================
# 5. Publication-level plotting 
#=======================================================

p <- ggplot(plot_df, aes(x = Age_num, y = target_norm)) +
  geom_jitter(width = 1.2, alpha = 0.4, size = 1, color = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, color = "#c0392b", fill = "#e74c3c", alpha = 0.2, linewidth = 0.8) +
  facet_wrap(~ Tissue, scales = "free_y", ncol = 5) + 
  geom_text(
    data = lm_results,
    aes(x = x_text, y = y_text, label = label),
    vjust = 1, hjust = 0, size = 2.8, color = "black", fontface = "italic",
    inherit.aes = FALSE
  ) +
  labs(
    x = "Age (Years)",
    y = expression(paste(italic("BHLHE41"), " Expression (log", " CPM)")),
    title = "Age-associated expression of BHLHE41 across human tissues"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#ecf0f1", color = "white"),
    strip.text = element_text(size = 9, face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text = element_text(color = "black")
  )

print(p)



library(ComplexHeatmap)
library(circlize)

#=======================================================
# 1. Prepare and organize data matrix
#=======================================================

# (Assume mat_z is the Tissue x Age_num matrix generated in previous steps, ordered by clustering)
# Name matrix columns to display as age groups on the heatmap
colnames(mat_z) <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")

# Ensure the order of statistical results (lm_results) exactly matches the row names of the matrix
lm_matched <- lm_results[match(rownames(mat_z), lm_results$Tissue), ]

# Build text labels displayed on the right: tissue name, Beta value, FDR value
tissue_names <- rownames(mat_z)
beta_text <- sprintf("%.3f", lm_matched$estimate)
# Concatenate FDR and significance stars, maintaining scientific notation and alignment
fdr_text <- paste0(sprintf("%.2e", lm_matched$FDR), " ", lm_matched$sig_star)

#=======================================================
# 2. Define complex statistical and tissue annotations on the right
#=======================================================

right_annot <- rowAnnotation(
  
  # 1. Tissue name column
  Tissue = anno_text(tissue_names, 
                     just = "left",          
                     location = 0,           
                     gp = gpar(fontsize = 10, col = "black")),
  
  # 2. Beta value column
  Beta = anno_text(beta_text, 
                   just = "right",          
                   location = 1,           
                   gp = gpar(fontsize = 10, col = "black")),
  
  # 3. FDR value column
  "P-value" = anno_text(fdr_text, 
                        just = "left",          
                        location = 0,           
                        gp = gpar(fontsize = 10, col = "black", fontface = "italic")),
  
  # Key settings: ensure titles are at the top, without rotation
  annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
  annotation_name_side = "top", 
  annotation_name_rot = 0,
  
  # Spacing between columns
  gap = unit(4, "mm") 
)

#=======================================================
# 3. Define new heatmap color mapping (up-regulation red, down-regulation green)
#=======================================================

# Z-score: negative values indicate relative down-regulation (green/blue), 0 indicates average level (white), positive values indicate relative up-regulation (red)
col_fun <- colorRamp2(
  breaks = c(-2, 0, 2), 
  colors = c("#00AFBB", "white", "#F8766D")
)

#=======================================================
# 4. Draw and output final clustered heatmap
#=======================================================

ht <- Heatmap(
  mat_z,
  name = "Row\nZ-score",             
  col = col_fun,                     # Apply new color scheme
  
  # Clustering settings
  cluster_rows = TRUE,               
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE,           
  
  # Dendrogram and row name settings
  row_dend_side = "left",            
  show_row_names = FALSE,            # Turn off default row names on heatmap, replaced by annotation column
  row_dend_width = unit(2, "cm"),    
  
  # Add complex annotation on the right
  right_annotation = right_annot,
  
  # Font detail settings
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 0,             
  
  # Title
  column_title = "Age-associated expression of BHLHE41 across tissues (GTEx v10)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # Cell borders
  rect_gp = gpar(col = "white", lwd = 1)
)

# 【Core fix】：padding order is c(top, right, bottom, left)
# Increase top padding from 2mm to 18mm to "push out" hidden annotation headers!
draw(ht, padding = unit(c(18, 10, 2, 2), "mm"))
