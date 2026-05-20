library(dplyr)
library(scales)
library(cmapR)
library(GEOquery)
library(edgeR)
library(broom)
library(limma)
library(ggplot2)
# 仅保留用于出版级绘图

# 假设数据已经加载：tissue_exp, sample_annot, subject_annot

expr_mat <- tissue_exp@mat    
row_ann  <- tissue_exp@rid      
col_ann  <- tissue_exp@cid

# 目标基因
target_gene_id <- "ENSG00000123095.6"
Target_gene <- "BHLHE41"

#=======================================================
# 2. 构建 metadata (Base R 语法)
#=======================================================

# 提取并重命名 sample_info
# 注意：强烈建议加入缺血时间 SMTSISCH 作为协变量
sample_info <- sample_annot[, c("SAMPID", "SMTSD", "SMRIN", "SMTSISCH")]
colnames(sample_info) <- c("SAMPID", "Tissue", "RIN", "IschemicTime")
sample_info$SUBJID <- substr(sample_info$SAMPID, 1, 10)

# 提取并重命名 subject_info
subject_info <- subject_annot[, c("SUBJID", "AGE", "SEX")]
colnames(subject_info) <- c("SUBJID", "Age", "Sex")

# 合并 metadata 并过滤列注释中存在的样本
meta_full <- merge(sample_info, subject_info, by = "SUBJID")
meta_full <- meta_full[meta_full$SAMPID %in% col_ann, ]

# 年龄映射
age_map <- c(
  "20-29" = 25, "30-39" = 35, "40-49" = 45,
  "50-59" = 55, "60-69" = 65, "70-79" = 75
)
meta_full$Age_num <- age_map[as.character(meta_full$Age)]

# 清理表达矩阵
expr_mat <- as.matrix(expr_mat)
expr_mat <- expr_mat[!is.na(rownames(expr_mat)), !is.na(colnames(expr_mat))]

#=======================================================
# 3. 逐组织分析（含 QC 与 协变量校正）
#=======================================================

results_list <- list()
plot_data_list <- list()

tissues <- unique(meta_full$Tissue)

for (tissue in tissues) {
  cat("Processing:", tissue, "\n")
  
  #----- QC -----
  is_brain <- grepl("^Brain-", tissue)
  
  # 提取该组织的样本，并剔除年龄或 RIN 缺失的值
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
  
  #----- 提取表达矩阵与过滤 -----
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
  
  #----- 提取目标基因的表达数据 -----
  target_voom <- v$E[target_gene_id, ]
  
  df_plot <- meta_t
  # 使用 match 确保样本对齐
  df_plot$target_norm <- target_voom[match(df_plot$SAMPID, colnames(v$E))]
  
  #----- 多元线性回归 (Top Journal Standard) -----
  # 动态判断性别是否可以作为协变量 (排除卵巢、前列腺等)
  has_sex_covariate <- length(unique(df_plot$Sex)) > 1
  
  if (has_sex_covariate) {
    fit <- lm(target_norm ~ Age_num + Sex + RIN, data = df_plot) 
  } else {
    fit <- lm(target_norm ~ Age_num + RIN, data = df_plot)
  }
  
  # 提取模型结果中的 Age_num 统计量
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
# 4. 合并结果与 FDR 校正 (Base R 实现)
#=======================================================

if (length(results_list) == 0) stop("No tissues passed QC.")

# 合并列表为数据框
lm_results <- do.call(rbind, results_list)
plot_df <- do.call(rbind, plot_data_list)

# 计算 FDR 并添加显著性星号
lm_results$FDR <- p.adjust(lm_results$p.value, method = "BH")
lm_results$sig_star <- ifelse(lm_results$FDR < 0.001, "***",
                              ifelse(lm_results$FDR < 0.01, "**",
                                     ifelse(lm_results$FDR < 0.05, "*", "")))

# 构建绘图文本标签
lm_results$label <- paste0("N = ", lm_results$N, 
                           "\nβ = ", signif(lm_results$estimate, 3), 
                           "\nFDR = ", signif(lm_results$FDR, 3), " ", lm_results$sig_star)

# 动态计算标签的 x, y 坐标以防止遮挡 (使用 Base R 的 aggregate)
y_pos <- aggregate(target_norm ~ Tissue, data = plot_df, FUN = function(x) max(x, na.rm = TRUE) * 1.05)
colnames(y_pos) <- "y_text"

x_pos <- aggregate(Age_num ~ Tissue, data = plot_df, FUN = function(x) min(x, na.rm = TRUE))
colnames(x_pos) <- "x_text"

# 将坐标合并到 lm_results
# 1. 计算 y 坐标并强制指定所有列名
y_pos <- aggregate(target_norm ~ Tissue, data = plot_df, FUN = function(x) max(x, na.rm = TRUE) * 1.05)
colnames(y_pos) <- c("Tissue", "y_text") # 强制确保第一列叫 Tissue，第二列叫 y_text

# 2. 计算 x 坐标并强制指定所有列名
x_pos <- aggregate(Age_num ~ Tissue, data = plot_df, FUN = function(x) min(x, na.rm = TRUE))
colnames(x_pos) <- c("Tissue", "x_text") # 强制确保第一列叫 Tissue，第二列叫 x_text

# 3. 合并 x 和 y 坐标
pos_df <- merge(x_pos, y_pos, by = "Tissue")

# 4. 将坐标合并回统计结果 (检查 lm_results 是否也有 Tissue 列，正常情况下是有的)
lm_results <- merge(lm_results, pos_df, by = "Tissue")

#=======================================================
# 5. 出版级绘图 
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
# 1. 准备和整理数据矩阵
#=======================================================

# (确保 mat_z 是之前步骤生成的 Tissue x Age_num 矩阵，按聚类顺序排序)
# 为矩阵列命名，使其在热图上显示为年龄段
colnames(mat_z) <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")

# 确保统计结果 (lm_results) 的顺序与矩阵的行名完全一一对应
lm_matched <- lm_results[match(rownames(mat_z), lm_results$Tissue), ]

# 构建右侧显示的文本标签：组织名称、Beta 值、FDR 值
tissue_names <- rownames(mat_z)
beta_text <- sprintf("%.3f", lm_matched$estimate)
# 拼接 FDR 和显著性星号，保持科学计数法和对齐
fdr_text <- paste0(sprintf("%.2e", lm_matched$FDR), " ", lm_matched$sig_star)

#=======================================================
# 2. 定义右侧的复杂统计学和组织注释
#=======================================================

right_annot <- rowAnnotation(
  
  # 1. 组织名称列
  Tissue = anno_text(tissue_names, 
                     just = "left",          
                     location = 0,           
                     gp = gpar(fontsize = 10, col = "black")),
  
  # 2. Beta 值列
  Beta = anno_text(beta_text, 
                   just = "right",          
                   location = 1,           
                   gp = gpar(fontsize = 10, col = "black")),
  
  # 3. FDR 值列
  "P-value" = anno_text(fdr_text, 
                        just = "left",          
                        location = 0,           
                        gp = gpar(fontsize = 10, col = "black", fontface = "italic")),
  
  # 关键设置：确保标题在顶部，且不旋转
  annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
  annotation_name_side = "top", 
  annotation_name_rot = 0,
  
  # 列之间的间距
  gap = unit(4, "mm") 
)

#=======================================================
# 3. 定义全新热图颜色映射 (上调红，下调绿)
#=======================================================

# Z-score: 负值表示相对下调(绿/蓝)，0表示平均水平(白)，正值表示相对上调(红)
col_fun <- colorRamp2(
  breaks = c(-2, 0, 2), 
  colors = c("#00AFBB", "white", "#F8766D")
)

#=======================================================
# 4. 绘制并输出终极版聚类热图
#=======================================================

ht <- Heatmap(
  mat_z,
  name = "Row\nZ-score",             
  col = col_fun,                     # 应用新的配色
  
  # 聚类设置
  cluster_rows = TRUE,               
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE,           
  
  # 树状图与行名设置
  row_dend_side = "left",            
  show_row_names = FALSE,            # 热图自带行名关闭，由注释列代替
  row_dend_width = unit(2, "cm"),    
  
  # 添加右侧复杂注释
  right_annotation = right_annot,
  
  # 字体细节设置
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 0,             
  
  # 标题
  column_title = "Age-associated expression of BHLHE41 across tissues (GTEx v10)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # 单元格边框
  rect_gp = gpar(col = "white", lwd = 1)
)

# 【核心修复】：padding 的顺序是 c(上, 右, 下, 左)
# 将顶部的 padding 从 2mm 增加到 18mm，把隐藏的注释抬头"逼"出来！
draw(ht, padding = unit(c(18, 10, 2, 2), "mm"))