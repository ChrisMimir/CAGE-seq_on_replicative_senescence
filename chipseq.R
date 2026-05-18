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
