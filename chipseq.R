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

# 1. 导入 "Gained" 位点 (Treatment > Control)
gained_sites_path <- "diff_sites_final_c3.0_cond1.bed"
gained_sites <- import(gained_sites_path)

# 2. 导入 "Lost" 位点 (Control > Treatment)
lost_sites_path <- "diff_sites_final_c3.0_cond2.bed"
lost_sites <- import(lost_sites_path)
print(gained_sites)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno_gained <- annotatePeak(gained_sites, tssRegion=c(-3000, 3000), TxDb=txdb)
#peakAnno_lost <- annotatePeak(lost_sites, tssRegion=c(-2000,2000), TxDb=txdb)

plotAnnoBar(peakAnno_gained)

# 1. 将 peakAnno_gained 对象转换为 data.frame，这是最容易操作的
peak_df <- as.data.frame(peakAnno_gained)

# 查看一下列名，你会看到 'annotation' 和 'geneId'
head(peak_df)
#-----------
# 2. **关键检查：** 确保 score 列是数值型
peak_df$score <- as.numeric(peak_df$score)

# 3. 绘制小提琴图和叠加箱线图
p_violin_all <- ggplot(peak_df, aes(x = 1, y = score)) +
  
  # 绘制小提琴图
  geom_violin(fill = "#8A2BE2",      # 紫蓝色填充
              alpha = 0.7, 
              width = 0.5) +
  
  # 叠加箱线图，显示中位数和四分位数
  geom_boxplot(width = 0.1, 
               fill = "white", 
               color = "black", 
               outlier.shape = NA) + # 隐藏离群点，让小提琴图展示离群点
  
  # 增加中位数线 (深红色点和线)
  stat_summary(fun = median, 
               geom = "point", 
               size = 3, 
               color = "darkred") +
  stat_summary(fun = median, 
               geom = "line", 
               size = 1, 
               color = "darkred",
               group = 1) +
  
  # 添加标题和标签
  labs(title = "Bindng score of Gained Sites of SP1 in IR HUVECs",
       y = "Score",
       x = "") +
  
  # 主题美化：移除 X 轴刻度和标签
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold") # 标题居中加粗
  )

# 显示图形
print(p_violin_all)

# 2. 计算 score 列的中位数
# 使用 na.rm = TRUE 确保忽略任何缺失值（NA），防止计算失败
median_score_all_sites <- median(peak_df$score, na.rm = TRUE)

# 3. 输出中位数
print(paste("所有 Gained Sites (peak_df) 中 Score 的中位数 (Median) 是:", median_score_all_sites))

# 4. (可选) 计算所有分位数 (Q1, Median, Q3)
quartiles_all_sites <- quantile(peak_df$score, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
print("所有 Gained Sites 的分位数 (Q1, Median, Q3):")
print(quartiles_all_sites)

#-----------
# 2. 筛选出所有 'annotation' 列为 "Promoter (<=1kb)" 的行
# (注意: 字符串必须完全匹配 'plotAnnoBar' 中显示的标签)

promoter_TSS_df <- peak_df
#promoter_TSS_df <- peak_df[peak_df$annotation =="Promoter (<=1kb)", ]

promoter_TSS_df <- peak_df[abs(peak_df$distanceToTSS)<= 1000, ]
#promoter_TSS_df <- peak_df[peak_df$score > 4, ]
# 3. 检查一下筛选结果
print(paste("总共 Gained sites:", nrow(peak_df)))
print(paste("其中 selected sites:", nrow(promoter_TSS_df)))
head(promoter_TSS_df)


# 3. (可选) 将这些TSS信息保存为 CSV 文件
write.csv(promoter_TSS_df, file = "Gained_gene_list.csv", row.names = FALSE)


# 1. 从筛选后的 data.frame ('promoter_1kb_df') 中提取 'geneId'
#    'geneId' 通常是 Entrez ID
promoter_entrez_ids <- promoter_TSS_df$geneId

# 2. 去除NA (有些 peak 可能没有注释到基因) 和 重复项
promoter_entrez_ids_unique <- unique(promoter_entrez_ids[!is.na(promoter_entrez_ids)])

print(paste("用于通路富集的唯一基因数量:", length(promoter_entrez_ids_unique)))

# 3. 运行 GO 富集分析
go_results_promoter <- enrichGO(gene          = promoter_entrez_ids_unique,
                                OrgDb         = org.Hs.eg.db,
                                keyType       = 'ENTREZID',
                                ont           = "BP", # BP, MF, CC
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff  = 0.05)

# 4. (可选) 运行 KEGG 富集分析
kegg_results_promoter <- enrichKEGG(gene         = promoter_entrez_ids_unique,
                                    organism     = 'hsa', # 'hsa' for human
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff  = 0.05)

# 5. 可视化富集结果

dotplot(go_results_promoter, showCategory=20)

dotplot(kegg_results_promoter, showCategory=20)

files <- list(
  Young   = "con_Y_peaks.narrowPeak",
  IR = "con_IR_peaks.narrowPeak"
)

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-3000, 3000), verbose=FALSE)

plotAnnoBar(peakAnnoList)

#--------------------------------------------------

#---------------------------with chipseq data-----------------------
# Ensure row names of res_up and s_sharp match
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

# 1. Three-way intersection: A & B & C
ABC <- length(intersect(intersect(A, B), C))

# 2. Pairwise intersections (only the intersection parts, excluding the three-way intersection):
# A & B (only) = (A ∩ B) - (A ∩ B ∩ C)
AB <- length(setdiff(intersect(A, B), C))
AC <- length(setdiff(intersect(A, C), B))
BC <- length(setdiff(intersect(B, C), A))

# 3. Unique parts (only themselves):
# A (only) = A - (A ∩ B) - (A ∩ C) + (A ∩ B ∩ C)
A_only <- length(setdiff(setdiff(A, B), C))
B_only <- length(setdiff(setdiff(B, A), C))
C_only <- length(setdiff(setdiff(C, A), B))

fit_data <- c(
  "ChIP-seq" = A_only,
  "SharpTSS_CAGE-seq" = B_only,
  "broadTSS_CAGE-seq" = C_only,
  
  "ChIP-seq&SharpTSS_CAGE-seq" = AB,
  "ChIP-seq&broadTSS_CAGE-seq" = AC,
  "SharpTSS_CAGE-seq&broadTSS_CAGE-seq" = BC,
  
  "ChIP-seq&SharpTSS_CAGE-seq&broadTSS_CAGE-seq" = ABC
)
#BiocManager::install("eulerr")
#library(eulerr)
# 2. Fit Euler model
fit <- euler(fit_data)
plot(fit,
     # Key settings:
     quantities = TRUE, # Display quantities (Count)
     
     # Color and style settings: 
     labels = set_names,
     # Fill colors: select three base colors
     fills = list(fill = c("#00B0F6", "#F8766D", "#53B400"), alpha = 0.6), 
     
     # Set the style for intersection color blending
     # eulerr blends colors in intersection areas by default; here we just customize the base colors
     # To show clear blended colors in intersection areas, the alpha value must be less than 1; we have set it to 0.6
     
     edges = list(col = "black", lex = 2),
     main = "Ovelapping of annotated genes"
)


# 1. ChIP-seq Gained Sites enrichment results
chip_go_results <- enrichGO(gene          = chip_entrez_ids_char,
                            OrgDb         = org.Hs.eg.db,
                            keyType       = 'ENTREZID',
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff  = 0.05)

# 2. S_CAGE-seq differentially expressed genes enrichment results
s_cage_go_results <- enrichGO(gene          = s_geneIDs_char,
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

s_kegg_chip <- enrichKEGG(
  gene          =  chip_entrez_ids_char,
  organism      = "hsa",       # KEGG abbreviation for the corresponding species, e.g., "hsa" for human
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

s_kegg_sharp <- enrichKEGG(
  gene          =  s_geneIDs_char,
  organism      = "hsa",       # KEGG abbreviation for the corresponding species, e.g., "hsa" for human
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

dotplot(s_kegg_chip, showCategory = 20, title = "KEGG Enrichment for SP1-binding genes") +
  theme_minimal()

dotplot(s_kegg_sharp, showCategory = 20, title = "KEGG Enrichment for sharp TSS derived genes") +
  theme_minimal()

# 1. Convert to data frame and filter pathways with Q < 0.05
chip_df <- as.data.frame(chip_go_results)
s_cage_df <- as.data.frame(s_cage_go_results)

# 2. Filter pathways with Q-value less than 0.05 (this step is optional if cutoff was set in enrichGO)
chip_sig_df <- chip_df[chip_df$qvalue < 0.05, ]
s_cage_sig_df <- s_cage_df[s_cage_df$qvalue < 0.05, ]
chip_sig_df 

# 1. Extract the set of significant pathway IDs for ChIP-seq
chip_pathway_ids <- chip_sig_df$ID

# 2. Extract the set of significant pathway IDs for S_CAGE-seq
s_cage_pathway_ids <- s_cage_sig_df$ID

print(paste("Total significant pathways for ChIP-seq:", length(chip_pathway_ids)))
print(paste("Total significant pathways for S_CAGE-seq:", length(s_cage_pathway_ids)))

# Define sets A and B
A <- chip_pathway_ids
B <- s_cage_pathway_ids
set_names <- c("ChIP-seq ", "SharpTSS_CAGE-seq ")

# 1. Calculate intersection: A & B
A_and_B <- length(intersect(A, B))

# 2. Unique parts:
A_only <- length(setdiff(A, B))
B_only <- length(setdiff(B, A))

# 3. Create the named vector required by eulerr
fit_data_pathways <- c(
  "ChIP-seq " = A_only,
  "SharpTSS_CAGE-seq " = B_only,
  "ChIP-seq &SharpTSS_CAGE-seq " = A_and_B
)

print("Pathway count statistics for Venn diagram:")
print(fit_data_pathways)

# 1. Fit Euler model
# Note: eulerr::euler must be used
fit_pathways <- eulerr::euler(fit_data_pathways)

# 2. Draw the plot
plot(fit_pathways,
     # Key settings: area proportional + display quantities
     quantities = TRUE, 
     
     # Labels and colors
     labels = set_names,
     fills = list(fill = c("#00C1AA", "#F8766D"), alpha = 0.6), # Use a different color combination
     edges = list(col = "black", lex = 2),
     
     # Adjust margins to prevent text overlap
     margin = 0.15,          
     labels.cex = 1.0,       # Set label text size
     quantities.cex = 1.2,   # Number size inside circles
     
     main = "Overlapping pathways (GO)"
)
# kegg enrichment

# 1. ChIP-seq Gained Sites KEGG enrichment results
chip_kegg_results <- enrichKEGG(gene         = chip_entrez_ids_char,
                                organism     = 'hsa', 
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05)

# 2. S_CAGE-seq differentially expressed genes KEGG enrichment results
s_cage_kegg_results <- enrichKEGG(gene         = s_geneIDs_char,
                                  organism     = 'hsa', 
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05)

# 1. Convert to data frame and filter pathways with Q < 0.05
chip_kegg_df <- as.data.frame(chip_kegg_results)
s_cage_kegg_df <- as.data.frame(s_cage_kegg_results)

# 2. Filter pathway IDs with Q-value less than 0.05
chip_kegg_ids <- chip_kegg_df[chip_kegg_df$qvalue < 0.05, ]$ID
s_cage_kegg_ids <- s_cage_kegg_df[s_cage_kegg_df$qvalue < 0.05, ]$ID

print(paste("Total significant KEGG pathways for ChIP-seq:", length(chip_kegg_ids)))
print(paste("Total significant KEGG pathways for S_CAGE-seq:", length(s_cage_kegg_ids)))

# 3. Define set names and calculate quantities
A <- chip_kegg_ids
B <- s_cage_kegg_ids
set_names <- c("ChIP-seq", "SharpTSS_CAGE-seq")

# A & B
A_and_B_kegg <- length(intersect(A, B))
# Unique to A
A_only_kegg <- length(setdiff(A, B))
# Unique to B
B_only_kegg <- length(setdiff(B, A))

# 4. Create the named vector required by eulerr
fit_data_kegg <- c(
  "ChIP-seq" = A_only_kegg,
  "SharpTSS_CAGE-seq" = B_only_kegg,
  "ChIP-seq&SharpTSS_CAGE-seq" = A_and_B_kegg
)

print("KEGG pathway count statistics for Venn diagram:")
print(fit_data_kegg)

# 1. Fit Euler model
# Note: explicitly specify function source using eulerr::euler
fit_kegg <- eulerr::euler(fit_data_kegg)

# 2. Draw the plot
plot(fit_kegg,
     # Key settings: area proportional + display quantities
     quantities = TRUE, 
     
     # Labels and colors
     labels = set_names,
     fills = list(fill = c("#00C1AA", "#F8766D"), alpha = 0.6), # Adopt KEGG-style green and red
     edges = list(col = "black", lex = 2),
     
     # Adjust margins to prevent text overlap
     margin = 0.15,          
     labels.cex = 1.0,       # Set label text size
     quantities.cex = 1.2,   # Number size inside circles
     
     main = "Overlapping pathways(KEGG)"
)
#------------------------------
# 1. Extract unique Entrez IDs for common genes (ensure no NAs)
ids_to_convert <- unique(common_sites_df$geneId[!is.na(common_sites_df$geneId)])
ids_to_convert <- as.character(ids_to_convert)
id_map <- mapIds(x = org.Hs.eg.db,
                 keys = ids_to_convert,
                 column = "SYMBOL",
                 keytype = "ENTREZID",
                 multiVals = "first") # If an Entrez ID corresponds to multiple Symbols, take only the first one

# 3. Convert the map into a data.frame
id_map_df <- data.frame(geneId = names(id_map), SYMBOL = id_map, stringsAsFactors = FALSE)

# 4. Merge Gene Symbols back into your common sites data frame
final_results_df <- merge(common_sites_df, id_map_df, by = "geneId", all.x = TRUE)

# 5. Sort descending by the 'score' column (assuming 'score' is your Fold Change or Log10 FDR)
# Note: Replace 'score' with the actual column name in your file (e.g., 'log2FC' or 'score')
# Assuming you want the highest Fold Change first; adjust sorting logic if you want the smallest FDR first.
final_results_sorted <- final_results_df[order(final_results_df$score, decreasing = TRUE), ]
print("Final results (filtered, converted to Gene Symbols, and sorted by score):")
# Display the first few rows, including the new SYMBOL column
head(final_results_sorted[, 
    c("SYMBOL", "geneId", "score", "annotation", "distanceToTSS")])

high_score <- final_results_sorted[(final_results_sorted$distanceToTSS > -10)&
                                     (final_results_sorted$distanceToTSS < 5),]
unique(high_score[, c("SYMBOL", "geneId")])
share_target <- unique(final_results_sorted[, c("SYMBOL", "geneId")])
share_target

set1_symbols <- as.character(gene_list$`1`)
set1_symbols2 <- as.character(share_target$SYMBOL)
common_symbols <- intersect(set1_symbols, set1_symbols2)
common_symbols

s_seqs <- s_sharp %>%
  rowRanges() %>%
  swapRanges() %>%
  promoters(upstream = 20, downstream = 10) %>%
  getSeq(bsg, .)

opts <- list(
  species = "9060",
  ID = "MA0636.1"        # ID = "MA0693.1" can also be used
)

motif_pfms <- getMatrixSet(JASPAR2020, opts)


motif_hits <- matchMotifs(
  motif_pfms,
  subject = s_seqs,
  out = "matches",     # Return match positions and P-values
  p.cutoff = 0.05
)

motif_logical <- motifMatches(motif_hits)  # sparse matrix
hit_idx <- which(motif_logical[, 1])
hit_tss <- s_sharp[hit_idx]
gene_hits <- rowData(hit_tss)$geneID

gene_symbols <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_hits,
  columns = "SYMBOL",
  keytype = "ENTREZID"     # Change to "ENTREZID" if geneID is Entrez
)

print(gene_symbols$SYMBOL)

SP1_genes <- setdiff(common_symbols, gene_symbols$SYMBOL)

SP1_genes
#----------------------------
BiocManager::install("ChIPseeker")
library(ChIPseeker)

# !!! Replace with your file path !!!
background_peak_path <- "diff_sites_final_c3.0_common.bed" 
common_peaks<- import(background_peak_path)

gained_sites_path <- "diff_sites_final_c3.0_cond1.bed"
gained_sites <- import(gained_sites_path)


all_chip_peaks<- c(gained_sites, common_peaks)
# Find common chromosome names
common_chroms <- intersect(seqlevels(gained_sites), seqlevels(all_cage_tss_gr))
common_chroms <- intersect(common_chroms, seqlevels(all_chip_peaks))

# Unify chromosome levels
seqlevels(gained_sites, pruning.mode="coarse") <- common_chroms
seqlevels(all_cage_tss_gr, pruning.mode="coarse") <- common_chroms
seqlevels(all_chip_peaks, pruning.mode="coarse") <- common_chroms


sharp_gr <- rowRanges(s_sharp)
broad_gr <- rowRanges(s_broad)

# 2b. Merge Sharp and Broad TSS coordinates
#all_cage_tss_gr <- c(sharp_gr, broad_gr)
all_cage_tss_gr <- broad_gr

# 2. Extend CAGE-seq TSSs (+/- 1kb)
cage_promoter_regions_gr <- resize(all_cage_tss_gr, width = 1000, fix = "center")


# 3. List B: TSS-associated Peaks (among all Peaks)
tss_overlaps_in_all <- findOverlaps(all_chip_peaks, cage_promoter_regions_gr)
tss_associated_peaks <- unique(queryHits(tss_overlaps_in_all))
m <- length(tss_associated_peaks) # K: Number of TSS-associated peaks in the universe

# 4. List A: Differentially bound Peaks (Gained Sites)
#    To ensure gained_sites are strictly contained in all_chip_peaks, we first take the intersection (safe approach)
gained_sites_in_universe <- intersect(gained_sites, all_chip_peaks)
gained_sites_count <- length(gained_sites_in_universe) # n: Sample size

# 5. Overlap: Peaks in Gained Sites that overlap with TSS (List X)
gained_tss_overlaps <- findOverlaps(gained_sites_in_universe, cage_promoter_regions_gr)
x <- length(unique(queryHits(gained_tss_overlaps))) # x: Overlap count

# Universe size
n_total_universe <- length(all_chip_peaks) 

# --- Construct matrix elements ---
# 1. Differentially bound AND TSS-associated (x)
overlap_gained_tss <- x

# 2. Differentially bound AND Non-TSS-associated (Unique to Gained Sites)
gained_non_tss <- gained_sites_count - overlap_gained_tss

# 3. Non-differentially bound AND TSS-associated (Background peaks associated with TSS, but not in Gained Sites)
# All TSS-associated peaks (m) - overlapping TSS-associated peaks (x)
background_tss_non_gained <- m - overlap_gained_tss

# 4. Neither (Background peaks - Total Gained Sites - Remaining TSS-associated peaks)
neither_gained_nor_tss <- n_total_universe - (overlap_gained_tss + gained_non_tss + background_tss_non_gained)

# Ensure all numbers are >= 0
if (any(c(overlap_gained_tss, gained_non_tss, background_tss_non_gained, neither_gained_nor_tss) < 0)) {
  stop("Negative numbers appeared in matrix calculation; please check if chromosome naming in GRanges objects is completely consistent.")
}


# 5. Construct the matrix
contingency_matrix_peak <- matrix(
  c(overlap_gained_tss, gained_non_tss,
    background_tss_non_gained, neither_gained_nor_tss),
  nrow = 2,
  dimnames = list(
    TSS_Association = c("TSS_Associated", "Non_TSS_Associated"),
    Diff_Peak_Type = c("Gained_Sites", "Non_Gained_Sites")
  )
)

# 6. Run Fisher's Exact Test
fisher_peak_result <- fisher.test(contingency_matrix_peak, alternative = "greater")


# --- 7. Output results ---
print("--- Peak Overlap (TSS +/- 1kb) Statistical Test ---")
print(paste("Percentage of TSS association in Gained Sites:", round(x / k * 100, 2), "%"))
print(paste("Percentage of TSS association in all Peaks:", round(m / n_total_universe * 100, 2), "%"))
print("\nContingency Matrix:")
print(contingency_matrix_peak)
print("\nFisher's Exact Test results:")
print(fisher_peak_result)
