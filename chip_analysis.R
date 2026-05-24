#----------load chip data--------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ChIPpeakAnno")
BiocManager::install("GenomicRanges")
BiocManager::install("DiffBind")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("ChIPseeker")
BiocManager::install("clusterProfiler")

library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GenomicRanges)
library(DiffBind)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(rtracklayer)
library(clusterProfiler)
library(ggplot2)

# 1. Import "Gained" sites (Treatment > Control)
gained_sites_path <- "diff_sites_final_c3.0_cond1.bed"
gained_sites <- import(gained_sites_path)

# 2. Import "Lost" sites (Control > Treatment)
lost_sites_path <- "diff_sites_final_c3.0_cond2.bed"
lost_sites <- import(lost_sites_path)
print(gained_sites)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno_gained <- annotatePeak(gained_sites, tssRegion=c(-3000, 3000), TxDb=txdb)
#peakAnno_lost <- annotatePeak(lost_sites, tssRegion=c(-2000,2000), TxDb=txdb)

plotAnnoBar(peakAnno_gained)

# 1. Convert peakAnno_gained object to data.frame, which is easiest to work with
peak_df <- as.data.frame(peakAnno_gained)

# Check column names, you will see 'annotation' and 'geneId'
head(peak_df)
#-----------
# 2. **Critical check:** Ensure score column is numeric
peak_df$score <- as.numeric(peak_df$score)

# 3. Create violin plot with overlaid boxplot
p_violin_all <- ggplot(peak_df, aes(x = 1, y = score)) +
  
  # Draw violin plot
  geom_violin(fill = "#8A2BE2",      # Purple-blue fill
              alpha = 0.7, 
              width = 0.5) +
  
  # Overlay boxplot to show median and quartiles
  geom_boxplot(width = 0.1, 
               fill = "white", 
               color = "black", 
               outlier.shape = NA) + # Hide outliers, let violin show them
  
  # Add median line (dark red points and line)
  stat_summary(fun = median, 
               geom = "point", 
               size = 3, 
               color = "darkred") +
  stat_summary(fun = median, 
               geom = "line", 
               size = 1, 
               color = "darkred",
               group = 1) +
  
  # Add title and labels
  labs(title = "Binding score of Gained Sites of SP1 in IR HUVECs",
       y = "Score",
       x = "") +
  
  # Theme beautification: remove X-axis ticks and labels
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold title
  )

# Display plot
print(p_violin_all)

# 2. Calculate median of score column
# Use na.rm = TRUE to ignore any missing values (NA), prevent calculation failure
median_score_all_sites <- median(peak_df$score, na.rm = TRUE)

# 3. Output median
print(paste("Median Score for all Gained Sites (peak_df) is:", median_score_all_sites))

# 4. (Optional) Calculate all quartiles (Q1, Median, Q3)
quartiles_all_sites <- quantile(peak_df$score, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
print("Quartiles for all Gained Sites (Q1, Median, Q3):")
print(quartiles_all_sites)

#-----------
# 2. Filter rows where 'annotation' column is "Promoter (<=1kb)"
# (Note: string must exactly match the label shown in 'plotAnnoBar')

promoter_TSS_df <- peak_df
#promoter_TSS_df <- peak_df[peak_df$annotation =="Promoter (<=1kb)", ]

#promoter_TSS_df <- peak_df[abs(peak_df$distanceToTSS)<= 1000, ]
#promoter_TSS_df <- peak_df[peak_df$score > 4, ]
# 3. Check filtering results
print(paste("Total Gained sites:", nrow(peak_df)))
print(paste("Selected sites:", nrow(promoter_TSS_df)))
head(promoter_TSS_df)


# 3. (Optional) Save TSS information as CSV file
write.csv(promoter_TSS_df, file = "Gained_gene_list.csv", row.names = FALSE)


#------------------ChIP data and CAGE data---------------------------------
# Ensure required packages are installed (if not already)
# BiocManager::install(c("org.Hs.eg.db", "ggplot2", "tidyr", "dplyr"))
library(org.Hs.eg.db) # Assuming human data
library(AnnotationDbi)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stats) # Fisher's Exact Test is in stats package
library(eulerr) 
library(cowplot) # Ensure loaded as_ggplot
library(ggplot2) # Ensure loaded annotate



# Ensure `res_up` and `s_sharp` row names can match
gained_genes_df <- read_csv("C:/Users/Jiaqi/Documents/Code/R/chipseq/Gained_gene_list.csv")
chip_entrez_ids <- unique(gained_genes_df$geneId[!is.na(gained_genes_df$geneId)])
common_entrez_ids <- intersect(chip_entrez_ids, s_geneIDs)
common_sites_df <- gained_genes_df[gained_genes_df$geneId %in% common_entrez_ids, ]
num_common <- length(common_entrez_ids)
chip_only_entrez_ids <- setdiff(chip_entrez_ids, s_geneIDs)
num_chip_only <- length(chip_only_entrez_ids)
cage_only_entrez_ids <- setdiff(s_geneIDs, chip_entrez_ids)
num_cage_only <- length(cage_only_entrez_ids)


chip_entrez_ids_char <- as.character(chip_entrez_ids)
s_geneIDs_char <- as.character(s_geneIDs)
b_geneIDs_char <- as.character(s_bgeneIDs)

# Define set names
set_names <- c("ChIP-seq", "SharpTSS_CAGE-seq", "broadTSS_CAGE-seq")

A <- chip_entrez_ids_char
B <- s_geneIDs_char
C <- b_geneIDs_char

# --- Define Universe (background gene set) ---
# Get all Entrez IDs from annotation package
universe_entrez_ids <- keys(org.Hs.eg.db, keytype="ENTREZID")

# Ensure all gene lists (A, B, C) are character type and only contain genes in Universe
A <- as.character(A)
B <- as.character(B)
C <- as.character(C)
U <- as.character(universe_entrez_ids)

# Restrict A, B, C to Universe range (this step is critical to avoid negative numbers in matrix)
A <- intersect(A, U)
B <- intersect(B, U)
C <- intersect(C, U)

# Ensure Universe contains all genes used for calculation
# Strictly speaking, Universe should be all genes detectable by these three methods.
# But using the entire species gene set (U) is the most standard approach.
U_size <- length(U)

# --- Helper function: Run Fisher's Test ---
run_fisher_test <- function(list1, list2, universe_size) {
  overlap <- length(intersect(list1, list2))
  only_list1 <- length(setdiff(list1, list2))
  only_list2 <- length(setdiff(list2, list1))
  
  # Neither = total Universe - (only list1 + only list2 + overlap)
  neither <- universe_size - (overlap + only_list1 + only_list2)
  
  # Build 2x2 matrix (note order when byrow=TRUE)
  matrix_data <- matrix(c(overlap, only_list1, only_list2, neither), 
                        nrow = 2, byrow = TRUE, 
                        dimnames = list(
                          In_A = c("Yes", "No"),
                          In_B = c("Yes", "No")
                        ))
  
  # Run test (greater tests for enrichment)
  test_result <- fisher.test(matrix_data, alternative = "greater")
  
  return(list(
    overlap = overlap,
    list1_unique = only_list1,
    p_value = test_result$p.value,
    odds_ratio = test_result$estimate,
    list1_size = length(list1),
    list2_size = length(list2)
  ))
}


# --- Run A vs B (ChIP vs Sharp) ---
result_AB <- run_fisher_test(A, B, U_size)
result_AB
# --- Run A vs C (ChIP vs Broad) ---
result_AC <- run_fisher_test(A, C, U_size)
result_AC

# --- Result summary (for plotting) ---
results_df <- data.frame(
  Comparison = c("ChIP-seq & SharpTSS", "ChIP-seq & BroadTSS"),
  Overlap_Count = c(result_AB$overlap, result_AC$overlap),
  ChIP_Unique = c(result_AB$list1_unique, result_AC$list1_unique),
  P_value = c(result_AB$p_value, result_AC$p_value),
  ChIP_Total = c(result_AB$list1_size, result_AC$list1_size)
)

# 1. Convert data to long format (required for Stacked Bar Plot)
plot_data <- results_df %>%
  select(Comparison, Overlap_Count, ChIP_Unique) %>%
  pivot_longer(cols = c(Overlap_Count, ChIP_Unique),
               names_to = "Category",
               values_to = "Count") %>%
  mutate(Category = factor(Category, 
                           levels = c("ChIP_Unique", "Overlap_Count"),
                           labels = c("ChIP-seq Unique", "Overlap with TSS")))

# 2. Create Stacked Bar Plot
p_overlap <- ggplot(plot_data, aes(x = Comparison, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  
  # Label overlap Count (white text on blue area)
  geom_text(data = filter(plot_data, Category == "Overlap with TSS"),
            aes(label = Count), 
            position = position_stack(vjust = 0.5), 
            color = "white", fontface = "bold", size = 4) +
  
  # Label P-value (critical fix: add inherit.aes = FALSE)
  geom_text(data = results_df, 
            aes(x = Comparison, y = ChIP_Total, 
                # Format P-value and add stars
                label = paste0("P = ", format.pval(P_value, digits = 2), 
                               ifelse(P_value < 0.001, "***", ifelse(P_value < 0.01, "**", ifelse(P_value < 0.05, "*", "")))
                )
            ),
            vjust = -0.5, size = 4, color = "darkred",
            inherit.aes = FALSE # <--- Prevents it from looking for Category mapping
  ) + 
  
  # Title, labels and theme
  scale_fill_manual(values = c("ChIP-seq Unique" = "gray", "Overlap with TSS" = "#00B0F6"),
                    name = "ChIP Target Gene Status") +
  labs(title = "Overlap Analysis of ChIP-seq Target Genes with Sharp/Broad TSS Target Genes",
       y = "Number of ChIP-seq Target Genes (Gene Count)",
       x = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    # Ensure Y-axis has enough space for P-value labels
    plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm") 
  )

# Display plot
print(p_overlap)

# -----------------------------------------------------------
# Figure 1: ChIP-seq vs. SharpTSS target gene overlap
# -----------------------------------------------------------

# 1. Calculate set counts
A_only_B <- length(setdiff(A, B))
B_only_A <- length(setdiff(B, A))
overlap_AB <- length(intersect(A, B))

fit_data_AB <- c(
  "ChIP-seq" = A_only_B,
  "SharpTSS" = B_only_A,
  "ChIP-seq&SharpTSS" = overlap_AB
)

# 2. Fit Euler model
fit_AB <- euler(fit_data_AB)

# 3. Create plot
plot_AB <- plot(fit_AB,
                quantities = TRUE,
                labels = c("ChIP-seq", "SharpTSS"),
                fills = list(fill = c("#53B400", "#F8766D"), alpha = 0.6), # Blue and red
                edges = list(col = "black", lex = 2),
                margin = 0.15,
                labels.cex = 1.0,
                quantities.cex = 1.2,
                main = "Gene overlap of ChIP-seq and sharp high TSSs"
)

# Print Figure 1
print(plot_AB)


# -----------------------------------------------------------
# Figure 2: ChIP-seq vs. BroadTSS target gene overlap
# -----------------------------------------------------------

# 1. Calculate set counts
A_only_C <- length(setdiff(A, C))
C_only_A <- length(setdiff(C, A))
overlap_AC <- length(intersect(A, C))

fit_data_AC <- c(
  "ChIP-seq" = A_only_C,
  "BroadTSS" = C_only_A,
  "ChIP-seq&BroadTSS" = overlap_AC
)

# 2. Fit Euler model
fit_AC <- euler(fit_data_AC)

# 3. Create plot
plot_AC <- plot(fit_AC,
                quantities = TRUE,
                labels = c("ChIP-seq", "BroadTSS"),
                fills = list(fill = c("#53B400", "#F8766D"), alpha = 0.6), # Blue and green
                edges = list(col = "black", lex = 2),
                margin = 0.15,
                labels.cex = 1.0,
                quantities.cex = 1.2,
                main = "Gene overlap of ChIP-seq and broad high TSSs"
)

# Print Figure 2
print(plot_AC)



#----------------------


# Ensure KEGGREST package is installed (if needed)
# BiocManager::install("KEGGREST") 
# library(KEGGREST) # or clusterProfiler internal method

# Run a lenient KEGG enrichment to get all available pathway IDs as background
# Using U (universe_entrez_ids) defined earlier
# Assuming U (universe_entrez_ids) is already defined and contains all Entrez IDs
all_pathways_results <- enrichKEGG(gene = U,
                                   organism = 'hsa',
                                   pvalueCutoff = 1,  # No P-value cutoff
                                   qvalueCutoff = 1) # No Q-value cutoff

# Extract all available KEGG pathway IDs as background set
universe_kegg_ids <- as.character(all_pathways_results$ID)
U_kegg <- unique(universe_kegg_ids)
U_kegg_size <- length(U_kegg)

print(paste("KEGG pathway background set (Universe) size:", U_kegg_size))


s_kegg_chip <- enrichKEGG(
  gene          =  A,
  organism      = "hsa",       # KEGG abbreviation for species, e.g., "hsa" for human
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

s_kegg_sharp <- enrichKEGG(
  gene          =  B,
  organism      = "hsa",       # KEGG abbreviation for species, e.g., "hsa" for human
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

dotplot(s_kegg_chip, showCategory = 20, title = "KEGG Enrichment for SP1-binding genes") +
  theme_minimal()

dotplot(s_kegg_sharp, showCategory = 20, title = "KEGG Enrichment for sharp TSS derived genes") +
  theme_minimal()

# 1. Convert to data frame and filter pathways with Q < 0.05
chip_df <- as.data.frame(s_kegg_chip)
s_cage_df <- as.data.frame(s_kegg_sharp)

# 2. Filter pathways with Q-value less than 0.05 (if you set cutoff in enrichGO, this step is optional)
chip_sig_df <- chip_df[chip_df$qvalue < 0.05, ]
s_cage_sig_df <- s_cage_df[s_cage_df$qvalue < 0.05, ]

# 1. Extract ChIP-seq significant pathway ID set
chip_pathway_ids <- chip_sig_df$ID

# 2. Extract S_CAGE-seq significant pathway ID set
s_cage_pathway_ids <- s_cage_sig_df$ID

print(paste("Total ChIP-seq significant pathways:", length(chip_pathway_ids)))
print(paste("Total S_CAGE-seq significant pathways:", length(s_cage_pathway_ids)))

# A: chip_kegg_ids (ChIP-seq significant pathways)
# B: s_cage_kegg_ids (SharpTSS significant pathways)

# 3. Define set names and calculate counts
A_kegg_chip <- chip_kegg_ids
B_kegg_sharp <- s_cage_kegg_ids

# Ensure A and B are within U_kegg range (safe practice)
A_in_U <- intersect(A_kegg_chip, U_kegg)
B_in_U <- intersect(B_kegg_sharp, U_kegg)

# 1. Calculate overlap (both A and B have)
overlap_kegg <- length(intersect(A_in_U, B_in_U)) # x

# 2. Calculate A unique (A has, B does not)
A_only_kegg <- length(setdiff(A_in_U, B_in_U)) # n1 - x

# 3. Calculate B unique (B has, A does not)
B_only_kegg <- length(setdiff(B_in_U, A_in_U)) # n2 - x

# 4. Calculate neither (in Universe, but neither A nor B have)
neither_kegg <- U_kegg_size - (overlap_kegg + A_only_kegg + B_only_kegg)

# 5. Build matrix
contingency_matrix_kegg <- matrix(
  c(overlap_kegg, A_only_kegg,
    B_only_kegg, neither_kegg),
  nrow = 2,
  dimnames = list(
    In_ChIP_Pathways = c("Yes", "No"),
    In_SharpTSS_Pathways = c("Yes", "No")
  )
)

# 6. Run Fisher's Exact Test
fisher_kegg_result <- fisher.test(contingency_matrix_kegg, alternative = "greater")


# --- 7. Output results ---
print("--- KEGG Pathway Overlap Fisher's Exact Test Results ---")
print("Contingency Matrix (pathway counts):")
print(contingency_matrix_kegg)
print("\nFisher's Exact Test Results:")
print(fisher_kegg_result)

#-----------------------------------------

# Reuse previous format_p_value function
format_p_value <- function(p_value) {
  if (is.null(p_value) || length(p_value) == 0 || is.na(p_value)) {
    return("P-value = N/A")
  }
  if (p_value < 0.01) {
    return(paste("P-value =", formatC(p_value, format = "e", digits = 2)))
  } else {
    return(paste("P-value =", round(p_value, 3)))
  }
}

# 1. Fit Euler model (reuse your provided fit_data_kegg)
fit_kegg <- eulerr::euler(fit_data_kegg)

# 1. Fit Euler model
fit_kegg <- eulerr::euler(fit_data_kegg)

# 2. Create Euler plot (base plot)
plot_kegg_base <- plot(fit_kegg,
                       quantities = TRUE,
                       labels = set_names,
                       fills = list(fill = c("#00C1AA", "#F8766D"), alpha = 0.6),
                       edges = list(col = "black", lex = 2),
                       margin = 0.1,
                       labels.cex = 1.0,
                       quantities.cex = 1.2,
                       main = NULL # Leave empty, add via ggplot labs
)

# 3. Convert to ggplot object and add annotation (fixed)
plot_kegg_final <- ggdraw(plot_kegg_base) +
  # Add main title
  labs(title = "Overlapping Pathways (KEGG)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  
  # Add P-value annotation
  annotate("text", 
           x = 0.05, y = 0.05, 
           label = format_p_value(fisher_kegg_result$p.value), 
           color = "darkred", size = 4.5, hjust = 0
  )

# Print final plot
print(plot_kegg_final)



#---------------------

all_pathways_results_go <- enrichGO(gene  = U,
                                    OrgDb         = org.Hs.eg.db,
                                    keyType       = 'ENTREZID',
                                    ont           = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 1,
                                    qvalueCutoff  = 1)

# Extract all available GO pathway IDs as background set
universe_go_ids <- as.character(all_pathways_results_go$ID)
U_go <- unique(universe_go_ids)
U_go_size <- length(U_go)

print(paste("GO pathway background set (Universe) size:", U_go_size))


# 1. ChIP-seq Gained Sites enrichment results
chip_go_results <- enrichGO(gene          = A,
                            OrgDb         = org.Hs.eg.db,
                            keyType       = 'ENTREZID',
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff  = 0.05)

# 2. S_CAGE-seq differentially expressed gene enrichment results
s_cage_go_results <- enrichGO(gene          = B,
                              OrgDb         = org.Hs.eg.db,
                              keyType       = 'ENTREZID',
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff  = 0.05)

dotplot(chip_go_results, showCategory = 20, title = "GO Enrichment for SP1-binding genes") +
  theme_minimal()

dotplot(s_cage_go_results, showCategory = 20, title = "GO Enrichment for sharp TSS derived gene") +
  theme_minimal()

# 1. Convert to data frame and filter pathways with Q < 0.05
chip_df_go <- as.data.frame(chip_go_results)
s_cage_df_go <- as.data.frame(s_cage_go_results)

# 2. Filter pathways with Q-value less than 0.05 (if you set cutoff in enrichGO, this step is optional)
chip_sig_df_go <- chip_df_go[chip_df_go$qvalue < 0.05, ]
s_cage_sig_df_go <- s_cage_df_go[s_cage_df_go$qvalue < 0.05, ]

# 1. Extract ChIP-seq significant pathway ID set
chip_pathway_ids_go <- chip_sig_df_go$ID

# 2. Extract S_CAGE-seq significant pathway ID set
s_cage_pathway_ids_go <- s_cage_sig_df_go$ID

print(paste("Total ChIP-seq significant pathways:", length(chip_pathway_ids_go)))
print(paste("Total S_CAGE-seq significant pathways:", length(s_cage_pathway_ids_go)))


A_chip_go <- chip_pathway_ids_go
B_sharp_go <- s_cage_pathway_ids_go

# Ensure A and B are within U_kegg range (safe practice)
A_in_U <- intersect(A_chip_go, U_go)
B_in_U <- intersect(B_sharp_go, U_go)

# 1. Calculate overlap (both A and B have)
overlap_go <- length(intersect(A_in_U, B_in_U)) # x

# 2. Calculate A unique (A has, B does not)
A_only_go <- length(setdiff(A_in_U, B_in_U)) # n1 - x

# 3. Calculate B unique (B has, A does not)
B_only_go <- length(setdiff(B_in_U, A_in_U)) # n2 - x

# 4. Calculate neither (in Universe, but neither A nor B have)
neither_go <- U_go_size - (overlap_go + A_only_go + B_only_go)

# 5. Build matrix
contingency_matrix_go <- matrix(
  c(overlap_go, A_only_go,
    B_only_go, neither_go),
  nrow = 2,
  dimnames = list(
    In_ChIP_Pathways = c("Yes", "No"),
    In_SharpTSS_Pathways = c("Yes", "No")
  )
)

# 6. Run Fisher's Exact Test
fisher_go_result <- fisher.test(contingency_matrix_go, alternative = "greater")


# --- 7. Output results ---
print("--- GO Pathway Overlap Fisher's Exact Test Results ---")
print("Contingency Matrix (pathway counts):")
print(contingency_matrix_go)
print("\nFisher's Exact Test Results:")
print(fisher_go_result)

#-----------------------------------------

# Reuse previous format_p_value function
format_p_value <- function(p_value) {
  if (is.null(p_value) || length(p_value) == 0 || is.na(p_value)) {
    return("P-value = N/A")
  }
  if (p_value < 0.01) {
    return(paste("P-value =", formatC(p_value, format = "e", digits = 2)))
  } else {
    return(paste("P-value =", round(p_value, 3)))
  }
}


# 1. Fit Euler model
fit_go <- eulerr::euler(fit_data_pathways)

# 2. Create Euler plot (base plot)
plot_go_base <- plot(fit_go,
                     quantities = TRUE,
                     labels = set_names,
                     fills = list(fill = c("#00C1AA", "#F8766D"), alpha = 0.6),
                     edges = list(col = "black", lex = 2),
                     margin = 0.1,
                     labels.cex = 1.0,
                     quantities.cex = 1.2,
                     main = NULL # Leave empty, add via ggplot labs
)

# 3. Convert to ggplot object and add annotation (fixed)
plot_go_final <- ggdraw(plot_kegg_base) +
  # Add main title
  labs(title = "Overlapping Pathways (GO)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  
  # Add P-value annotation (use relative coordinates x=0.05, y=0.05 to place at bottom-left)
  annotate("text", 
           x = 0.05, y = 0.05, # Place at bottom-left corner of canvas
           label = format_p_value(fisher_go_result$p.value), 
           color = "darkred", size = 4.5, hjust = 0
  )

# Print final plot
print(plot_go_final)



overlapped_go_ids <- intersect(chip_go_ids, s_cage_go_ids)

if (length(overlapped_go_ids) == 0) {
  stop("Error: No overlapping significant GO pathways (Q < 0.05). Cannot create Dotplot for overlapping terms.")
}

# --- 2. Filter data frames to keep only overlapping pathways ---
# Filter ChIP-seq data
chip_overlap_df <- chip_go_df %>%
  filter(ID %in% overlapped_go_ids) %>%
  mutate(Group = "ChIP-seq") %>%
  # Convert GeneRatio to numeric form for plotting (e.g., "10/100" -> 0.1)
  mutate(GeneRatio_val = eval(parse(text = GeneRatio))) 

# Filter Sharp TSS data
s_cage_overlap_df <- s_cage_go_df %>%
  filter(ID %in% overlapped_go_ids) %>%
  mutate(Group = "SharpTSS_CAGE-seq") %>%
  mutate(GeneRatio_val = eval(parse(text = GeneRatio))) 

# --- 3. Combine data frames ---
combined_go_overlap_df <- bind_rows(chip_overlap_df, s_cage_overlap_df)

# --- 4. Sorting and factorization ---
# Sort by pathway name to ensure points are side by side in the plot
combined_go_overlap_df <- combined_go_overlap_df %>%
  # Sort by ChIP-seq qvalue, putting significant ones at the front
  arrange(match(ID, chip_overlap_df$ID), qvalue) %>%
  # Convert pathway name (Description) to factor to preserve plotting order
  mutate(Description = factor(Description, 
                              levels = unique(Description)))




#--------------------

# --- 1. Identify overlapping pathway IDs ---
overlapped_pathway_ids <- intersect(chip_sig_df$ID, s_cage_sig_df$ID)

if (length(overlapped_pathway_ids) == 0) {
  stop("Error: No overlapping significant pathways between ChIP-seq and Sharp TSS (Q < 0.05).")
}

print(paste("Total number of overlapping significant pathways:", length(overlapped_pathway_ids)))

# --- 2. Filter original data frames to keep only overlapping pathways ---
chip_overlap_df <- chip_sig_df %>%
  filter(ID %in% overlapped_pathway_ids) %>%
  # Add a Group label
  mutate(Group = "ChIP-seq")

s_cage_overlap_df <- s_cage_sig_df %>%
  filter(ID %in% overlapped_pathway_ids) %>%
  mutate(Group = "SharpTSS_CAGE-seq")

# --- 3. Identify Top 15 pathway IDs ---
# Sorting strategy: Sort based on ChIP-seq P.adjust value (can be adjusted to SharpTSS or average as needed)
top_n <- 20

# a) Sort by P.adjust (smaller P.adjust comes first)
chip_sorted_overlap <- chip_overlap_df %>%
  arrange(p.adjust)

# b) Extract Top N pathway IDs (ensure uniqueness)
top_overlap_ids <- head(chip_sorted_overlap$ID, top_n)

print(paste("Number of Top", top_n, "overlapping pathways to display:", length(top_overlap_ids)))

# --- 4. Final Top N dataset ---
final_chip_top_df <- chip_overlap_df %>%
  filter(ID %in% top_overlap_ids)

final_s_cage_top_df <- s_cage_overlap_df %>%
  filter(ID %in% top_overlap_ids)

# --- 5. Prepare plotting data: Convert Description to factor to maintain Top 15 order ---

# Get final plotting order (based on ChIP-seq sorting)
plot_order <- final_chip_top_df %>%
  arrange(p.adjust) %>%
  pull(Description)

# Set pathway name factor levels for both data frames
final_chip_top_df <- final_chip_top_df %>%
  mutate(Description = factor(Description, levels = rev(plot_order))) # rev() puts most significant at top

final_s_cage_top_df <- final_s_cage_top_df %>%
  mutate(Description = factor(Description, levels = rev(plot_order)))


plot_chip <- ggplot(final_chip_top_df, 
                    aes(x = GeneRatio, y = Description, color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradientn(colors = c("#F8766D", "#C49A00", "#53B400"), name = "Adj. P-value") +
  scale_size_continuous(range = c(3, 8), name = "Gene Count") +
  labs(title = paste0("ChIP-seq (Top ", length(plot_order), " Overlap GO-BP)")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    # --- Critical modification: Rotate X-axis text ---
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate 45 degrees, right-align text to tick marks
    # If 90-degree rotation is needed: angle = 90, hjust = 1, vjust = 0.5
  )

print(plot_chip)

---
  
  # --- Plot 2: Sharp TSS Top 15 overlapping pathways ---
  plot_s_cage <- ggplot(final_s_cage_top_df, 
                        aes(x = GeneRatio, y = Description, color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradientn(colors = c("#F8766D", "#C49A00", "#53B400"), name = "Adj. P-value") +
  scale_size_continuous(range = c(3, 8), name = "Gene Count") +
  labs(title = paste0("Sharp TSS (Top ", length(plot_order), " Overlap GO-BP)")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    # --- Critical modification: Rotate X-axis text ---
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(plot_s_cage)