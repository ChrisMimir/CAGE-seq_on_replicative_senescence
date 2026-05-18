#Extended Data Figure 1a
# Step 1: 提取 TPM 并整理
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

# Step 2: 计算平均表达（按 txType 和 Class）
tpm_summary <- tpm_long %>%
  group_by(txType, Class) %>%
  summarise(mean_logTPM = mean(logTPM, na.rm = TRUE), .groups = "drop")

# 可选：确保 Class 顺序正确
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
           # 1. 修正点：去除黑色边框
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
    # 2. 修正点：去除小图标题的背景/边框
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

# 提取归一化计数矩阵
expr_matrix <- assay(TERSE, "size_factor_normalized_counts")  # 替换为合适的 assay 名称

# 提取分组信息
group_info <- colData(TERSE)$Class  # 分组信息
group_info <- factor(group_info)
average_exp <- apply(expr_matrix, 1, function(x) tapply(x, group_info, mean, na.rm = TRUE))

# 转换平均表达矩阵为数据框
average_exp_df <- as.data.frame(t(average_exp)) 
# 给列名添加分组名称
colnames(average_exp_df) <- levels(group_info)
# 将平均表达量添加到 rowData
rowData(TERSE)$average_exp <- average_exp_df

rowData(TERSE)$type <- "none"
rowData(TERSE)$transposon <- "none"
overlaps <- findOverlaps(rowRanges(RSE), TE_DE)

# 提取重叠的行索引
TERSE_hits <- queryHits(overlaps)
TE_DE_hits <- subjectHits(overlaps)

# 合并匹配的注释（处理重复匹配）
aggregated_type <- tapply(TE_DE$type[TE_DE_hits], TERSE_hits, paste, collapse = ",")
aggregated_transposon <- tapply(TE_DE$transposon[TE_DE_hits], TERSE_hits, paste, collapse = ",")

# 将注释信息对齐到 TERSE 的行
aligned_type <- rep("none", nrow(TERSE))
aligned_transposon <- rep("none", nrow(TERSE))

# 按索引更新注释信息
aligned_type[as.numeric(names(aggregated_type))] <- aggregated_type
aligned_transposon[as.numeric(names(aggregated_transposon))] <- aggregated_transposon

# 将对齐的注释信息添加到 rowData
rowData(TERSE)$type <- aligned_type
rowData(TERSE)$transposon <- aligned_transposon

# 检查结果
head(rowData(TERSE))

TERSE <- TERSE[rowData(TERSE)$type != "none", ]

# 提取归一化计数矩阵
expr_matrix <- assay(TERSE, "size_factor_normalized_counts")  # 替换为您的归一化计数矩阵

# 提取分组信息
group_info <- colData(TERSE)$Class  # 替换为您的分组列名

# 将分组信息转化为因子
group_info <- factor(group_info)

# 计算每组表达总计数
total_expression <- colSums(expr_matrix)

# 按分组统计总计数
group_totals <- tapply(total_expression, group_info, sum)

# 转换为数据框格式
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
# 使用 varianceStabilizingTransformation
te_blind <- varianceStabilizingTransformation(te_blind, blind = TRUE)

plotPCA(te_blind, "Class")

TEA_seqs <- TEA %>% 
  rowRanges() %>%
  swapRanges() %>%
  unstrand() %>%
  add(500) %>% 
  getSeq(bsg, names = .)
transposon_data <- rowData(TEA)$transposon

# 用 strsplit 拆分，并填充缺失值为 NA
transposon_split <- lapply(strsplit(transposon_data, "\\|"), function(x) {
  if (length(x) < 3) {
    c(x, rep(NA, 3 - length(x)))  # 如果长度小于 3，补充 NA
  } else {
    x[1:3]  # 如果长度大于或等于 3，仅取前 3 列
  }
})

transposon_split <- do.call(rbind, transposon_split)
# 将拆分结果命名为三列
colnames(transposon_split) <- c("class", "family", "subfamily")

# 将矩阵直接存储在 transposon 字段中
rowData(TEA)$transposon <- transposon_split
rowData(TEA)$transposon

# 获取所有转座子类别
transposon_classes <- unique(rowData(TEA)$transposon[, "class"])

calculate_transposon_matrix <- function(expr_matrix, group_info, transposon_classes, transposon_annotation) {
  # 初始化结果矩阵
  result_matrix <- matrix(0, nrow = length(transposon_classes), ncol = length(levels(group_info)))
  rownames(result_matrix) <- transposon_classes
  colnames(result_matrix) <- levels(group_info)
  
  # 按组计算每个类别的加权表达量
  for (group in levels(group_info)) {
    # 筛选当前组的样本
    group_samples <- which(group_info == group)
    
    # 当前组的表达矩阵
    group_expr <- expr_matrix[, group_samples]
    
    # 计算每个转座子类别的总表达量
    for (class in transposon_classes) {
      # 筛选当前类别的 TSS
      class_tss <- which(transposon_annotation == class)
      
      # 计算加权表达量总和
      result_matrix[class, group] <- sum(group_expr[class_tss, ], na.rm = TRUE)
    }
  }
  
  return(result_matrix)
}

# 提取必要数据
expr_matrix <- assay(TEA, "size_factor_normalized_counts")  # 替换为归一化计数矩阵
group_info <- factor(colData(TEA)$Class)  # 替换为分组信息
transposon_annotation <- rowData(TEA)$transposon[, "class"]  # 转座子类别注释

# 计算结果矩阵
transposon_matrix <- calculate_transposon_matrix(expr_matrix, group_info, transposon_classes, transposon_annotation)

head(transposon_matrix)

transposon_df <- as.data.frame(transposon_matrix)
transposon_df <- transposon_df %>%
  rownames_to_column(var = "class") %>%
  pivot_longer(cols = -class, names_to = "Group", values_to = "Expression")

# 绘制堆叠柱状图
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
# 取 TPM 矩阵（或 normalized counts）
tpm_matrix <- assay(TERSE, "TPM")  # 或 size_factor_normalized_counts

# 提取样本信息
sample_info <- as.data.frame(colData(TERSE)) %>%
  rownames_to_column("sample_name") %>%
  select(sample_name, Class)

# 提取注释信息（转座子与位置类型）
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

# 确保 Class 顺序一致
te_summary$Class <- factor(te_summary$Class, 
                           levels = c("p03", "p09", "p12", "p14", "p16", "p18"))

# PCA 图一致配色
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
# 假设 RSE 含有 TSS 层面的 reads
RSE_tss <- subset(RSE, clusterType == "TSS")

# 设计矩阵：时间组（例如 6 个传代）
colData(RSE_tss)$group <- factor(colData(RSE_tss)$Class,
                                 levels = c("p03","p09","p12","p14","p16","p18"))

dds <- DESeqDataSet(RSE_tss, design = ~ group)

# 过滤低表达 TSS
dds <- dds[rowSums(counts(dds)) > 10, ]

# DESeq 归一化与差异分析
dds <- DESeq(dds)

# 示例：p03 vs p18
res <- results(dds, contrast = c("group", "p18", "p03"))
res <- lfcShrink(dds, contrast = c("group", "p18", "p03"), res = res, type = "ashr")

# 加上 symbol 注释
res$symbol <- rowData(RSE_tss)$symbol[match(rownames(res), rownames(RSE_tss))]

# 添加显著标记
res$DE <- ifelse(res$padj < 0.05 & res$log2FoldChange > 1, "Up",
                 ifelse(res$padj < 0.05 & res$log2FoldChange < -1, "Down", "NS"))

#提取所有比较结果
comparisons <- combn(levels(colData(RSE_tss)$group), 2, simplify = FALSE)

# 汇总每组差异结果
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
  
  res_df <- as.data.frame(res)   # ✅ 转为 data.frame
  res_df$comparison <- paste(comp[2], "vs", comp[1])
  res_df$gene_id <- rownames(res_df)
  res_df$DE <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                      ifelse(res_df$log2FoldChange > 0, "Up", "Down"), "NS")
  return(res_df)
})

dt <- do.call(rbind, summary_list)

# 全局统计
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

# 提取并排序
enh_plot_df <- enh_summary %>%
  filter(comparison %in% selected_comparisons) %>%
  mutate(comparison = factor(comparison, levels = selected_comparisons))

# 自定义颜色
my_colors <- c(
  "#F8766D",
  "#C49A00",
  "#53B400",
  "#00C1AA",
  "#00B0F6"
)

# 绘制
p <- ggplot(enh_plot_df, aes(x = comparison, y = enh_count, fill = comparison)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = my_colors) +
  labs(
    x = "Contrast",
    y = "Number of enhancers linked to DE TSSs"
  ) +
  theme_bw(base_size = 14) +
  theme(
    # 去除右侧与上方边框
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    # 去除网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # X轴文字倾斜
    axis.text.x = element_text(angle = 45, hjust = 1),
    # 去掉图例
    legend.position = "none"
  )

# 展示
print(p)
