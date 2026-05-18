if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

#删除所有包
# remove.packages(row.names(installed.packages()))

# install required packages
BiocManager::install("CAGEfightR")
BiocManager::install(c("pheatmap","viridis","magrittr","ggforce","ggthemes","tidyverse"))
BiocManager::install(c("GenomicRanges","SummarizedExperiment","GenomicFeatures","BiocParallel","InteractionSet","Gviz"))
BiocManager::install(c("DESeq2","limma","edgeR","statmod","BiasedUrn","sva"))
BiocManager::install("clusterProfiler")
BiocManager::install("ggseqlogo")
BiocManager::install(c("TFBSTools","motifmatchr","pathview"))
BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38","TxDb.Hsapiens.UCSC.hg38.knownGene","org.Hs.eg.db"))
BiocManager::install("JASPAR2020")
BiocManager::install("Mfuzz")
BiocManager::install("msigdbr")
BiocManager::install("fgsea")
BiocManager::install("GSEABase")
install.packages("data.table")

# CRAN packages for data manipulation and plotting
library(pheatmap)
library(viridis)
library(magrittr)
library(ggforce)
library(ggthemes)
library(tidyverse)
library(ggseqlogo)

# CAGEfightR and related packages
library(CAGEfightR)
library(GenomicRanges)
library(SummarizedExperiment)
library(GenomicFeatures)
library(BiocParallel)
library(InteractionSet)
library(Gviz)

# Bioconductor packages for differential expression
library(DESeq2)
library(limma)
library(edgeR)
library(statmod)
library(BiasedUrn)
library(sva)

# Bioconductor packages for enrichment analyses
library(TFBSTools)
library(motifmatchr)
library(pathview)
library(clusterProfiler)

# Bioconductor data packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(JASPAR2020)
library(Mfuzz)

#library(data.table)
# Rename these for easier access
bsg <- BSgenome.Hsapiens.UCSC.hg38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
odb <- org.Hs.eg.db

# Script-wide settings
theme_set(theme_light()) # White theme for ggplot2 figures
register(SnowParam(6)) # Parallel execution when possible


#append data to data_matrix
data_matrix <- data.frame()
data_matrix <- rbind (data_matrix, data.frame(" " = "s01",Class = "p03", name = "s01",BigWigPlus = "p3_plus1.bw", BigWigMinus = "p3_minus1.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s02",Class = "p03", name = "s02",BigWigPlus = "p3_plus2.bw", BigWigMinus = "p3_minus2.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s03",Class = "p03", name = "s03",BigWigPlus = "p3_plus3.bw", BigWigMinus = "p3_minus3.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s07",Class = "p09", name = "s07",BigWigPlus = "p9_plus7.bw", BigWigMinus = "p9_minus7.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s08",Class = "p09", name = "s08",BigWigPlus = "p9_plus8.bw", BigWigMinus = "p9_minus8.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s09",Class = "p09", name = "s09",BigWigPlus = "p9_plus9.bw", BigWigMinus = "p9_minus9.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s10",Class = "p12", name = "s10",BigWigPlus = "p12_plus10.bw", BigWigMinus = "p12_minus10.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s11",Class = "p12", name = "s11",BigWigPlus = "p12_plus11.bw", BigWigMinus = "p12_minus11.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s12",Class = "p12", name = "s12",BigWigPlus = "p12_plus12.bw", BigWigMinus = "p12_minus12.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s13",Class = "p14", name = "s13",BigWigPlus = "p14_plus13.bw", BigWigMinus = "p14_minus13.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s14",Class = "p14", name = "s14",BigWigPlus = "p14_plus14.bw", BigWigMinus = "p14_minus14.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s15",Class = "p14", name = "s15",BigWigPlus = "p14_plus15.bw", BigWigMinus = "p14_minus15.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s16",Class = "p16", name = "s16",BigWigPlus = "p16_plus16.bw", BigWigMinus = "p16_minus16.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s17",Class = "p16", name = "s17",BigWigPlus = "p16_plus17.bw", BigWigMinus = "p16_minus17.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s18",Class = "p16", name = "s18",BigWigPlus = "p16_plus18.bw", BigWigMinus = "p16_minus18.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s19",Class = "p18", name = "s19",BigWigPlus = "p18_plus19.bw", BigWigMinus = "p18_minus19.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s20",Class = "p18", name = "s20",BigWigPlus = "p18_plus20.bw", BigWigMinus = "p18_minus20.bw"))
data_matrix <- rbind (data_matrix, data.frame(" " = "s21",Class = "p18", name = "s21",BigWigPlus = "p18_plus21.bw", BigWigMinus = "p18_minus21.bw"))

#colnames(data.matrix)[1] <- ""
#rownames(data_matrix) <- NULL

rownames(data_matrix)<-data_matrix[ ,1]
data.matrix<-data_matrix[,-1]
data_matrix <- subset(data_matrix, select = -c(X.))

knitr::kable(data_matrix, 
             caption = "The initial design matrix for the nanotubes experiment")

# Importing CTSSs

bw_plus <- BigWigFileList(c(
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p3_plus1.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p3_plus2.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p3_plus3.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p9_plus7.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p9_plus8.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p9_plus9.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p12_plus10.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p12_plus11.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p12_plus12.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p14_plus13.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p14_plus14.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p14_plus15.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p16_plus16.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p16_plus17.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p16_plus18.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p18_plus19.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p18_plus20.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p18_plus21.bw"
))
bw_minus <- BigWigFileList(c(
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p3_minus1.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p3_minus2.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p3_minus3.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p9_minus7.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p9_minus8.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p9_minus9.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p12_minus10.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p12_minus11.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p12_minus12.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p14_minus13.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p14_minus14.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p14_minus15.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p16_minus16.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p16_minus17.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p16_minus18.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p18_minus19.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p18_minus20.bw",
  "C:/Users/Jiaqi/Documents/Code/R/aging/raw_data/p18_minus21.bw"
))


names(bw_plus) <- data.matrix$name
names(bw_minus) <- data.matrix$name

CTSSs <- quantifyCTSSs(plusStrand = bw_plus,
                       minusStrand = bw_minus,
                       genome = seqinfo(bsg),
                       design = data_matrix)
CTSSs # Get a summary
rowRanges(CTSSs) # Extract CTSS positions
assay(CTSSs, "counts") %>%   # Extract CTSS counts
  head

CTSSs <- CTSSs %>%
  calcTPM() %>%
  calcPooled()

# CTSSs <- calcSupport(CTSSs,
#                      inputAssay="counts",
#                      outputColumn="support",
#                      unexpressed=0)
# table(rowRanges(CTSSs)$support)
# CTSSs <- subset(CTSSs, support > 1)
# CTSSs <- calcTPM(CTSSs, totalTags="totalTags")
# CTSSs <- calcPooled(CTSSs)

#-----------------------QC------------------
calcTotalTPM <- function(x) rowSums(assay(x, "TPM"))
rowData(CTSSs)$total_tpm <- calcTotalTPM(CTSSs)
ctss_data <- as.data.frame(rowData(CTSSs))
freq_table <- table(ctss_data$total_tpm)
freq_df <- data.frame(
  total_tpm = as.numeric(names(freq_table)),
  count = as.numeric(freq_table)
)

p1 <- ggplot(freq_df %>% filter(total_tpm <= 5), 
             aes(x = total_tpm, y = count)) +
  geom_col(fill = "steelblue", alpha = 0.7, width = 0.1) +  
  ggtitle("Distribution of Total TPM per CTSS (Low Range)") +
  xlab("Total TPM (Pooled across all samples)") +
  ylab("Number of CTSS") +
  theme_minimal()

thresholds <- seq(0, 2, by = 0.1)  
n_ctss <- sapply(thresholds, function(t) {
  sum(ctss_data$total_tpm > t)
})

p2 <- ggplot(data.frame(Threshold = thresholds, N_CTSS = n_ctss), 
             aes(x = Threshold, y = N_CTSS)) +
  geom_line(color = "red", size = 1) +
  geom_point(color = "red", size = 2) +
  ggtitle("Number of CTSS remaining after TPM thresholding") +
  ylab("Number of CTSS") +
  theme_minimal()

print(p1)
print(p2)

result_df <- data.frame(Threshold = thresholds, N_CTSS = n_ctss)
result_df$Percent_Remaining <- round(result_df$N_CTSS / n_ctss[1] * 100, 1)

cat("CTSS ratio under different threshold:\n")
print(result_df)

result_df$Reduction <- c(0, -diff(result_df$N_CTSS))
cat("\n CTSS reduction:\n")
print(result_df[, c("Threshold", "N_CTSS", "Reduction")])

#------------------------------------------
simple_TCs <- clusterUnidirectionally(CTSSs, 
                                      pooledCutoff=0.2, 
                                      mergeDist=20)

#TCs <- quickTSSs(CTSSs)
TCs <- quantifyClusters(CTSSs,
                        clusters=simple_TCs,
                        inputAssay="counts")

TSSs <- TCs %>%
  calcTPM() %>%
  subsetBySupport(inputAssay = "TPM", 
                  unexpressed = 1, 
                  minSamples = 2)

#BCs <- quickEnhancers(CTSSs)
BCs <- clusterBidirectionally(CTSSs, balanceThreshold=0.95)
# Calculate number of bidirectional samples
BCs <- calcBidirectionality(BCs, samples=CTSSs)
# Summarize
table(BCs$bidirectionality)

BCs <- subset(BCs, bidirectionality > 1)

Bidirectional <- quantifyClusters(CTSSs,
                                  clusters=BCs,
                                  inputAssay="counts")
BCs <- subsetBySupport(Bidirectional, 
                       inputAssay = "counts", 
                       unexpressed = 0, 
                       minSamples = 2)

# Annotate with transcript IDs
TSSs <- assignTxID(TSSs, txModels = txdb, swap = "thick")
# Annotate with transcript context
TSSs <- assignTxType(TSSs, txModels = txdb, swap = "thick")

# Annotate with transcript context
BCs <- assignTxType(BCs, txModels = txdb, swap = "thick")
# Keep only non-exonic BCs as enhancer candidates
Enhancers <- subset(BCs, txType %in% c("intergenic", "intron"))



# Merging into a single dataset
# Clean colData
TSSs$totalTags <- NULL
Enhancers$totalTags <- NULL

# Clean rowData
rowData(TSSs)$balance <- NA
rowData(TSSs)$bidirectionality <- NA
rowData(Enhancers)$txID <- NA

# Add labels for making later retrieval easy
rowData(TSSs)$clusterType <- "TSS"
rowData(Enhancers)$clusterType <- "Enhancer"

RSE <- combineClusters(object1 = TSSs, 
                       object2 = Enhancers, 
                       removeIfOverlapping = "object1")

RSE <- calcTPM(RSE)

# Genomic analysis of TSSs and enhancers

# Genome track
axis_track <- GenomeAxisTrack()

# Annotation track
tx_track <- GeneRegionTrack(txdb, 
                            name = "Gene Models", 
                            col = NA,
                            fill = "bisque4", 
                            shape = "arrow", 
                            showId = TRUE)

# Extract 100 bp around the first TSS
plot_region <- RSE %>% 
  rowRanges() %>% 
  subset(clusterType == "TSS") %>% 
  .[1] %>%
  add(100) %>%
  unstrand()

# CTSS track
ctss_track <- CTSSs %>%
  rowRanges() %>%
  subsetByOverlaps(plot_region) %>%
  trackCTSS(name = "CTSSs")

# Cluster track
cluster_track <- RSE %>%
  subsetByOverlaps(plot_region) %>%
  trackClusters(name = "Clusters", 
                col = NA, 
                showId = TRUE)

# Plot tracks together
plotTracks(list(axis_track, 
                ctss_track,
                cluster_track,
                tx_track),
           from = start(plot_region), 
           to = end(plot_region), 
           chromosome = as.character(seqnames(plot_region)))

# Extract 100 bp around the first enhancer
plot_region <- RSE %>% 
  rowRanges() %>% 
  subset(clusterType == "Enhancer") %>% 
  .[1] %>%
  add(100) %>%
  unstrand()

# CTSS track
ctss_track <- CTSSs %>%
  rowRanges() %>%
  subsetByOverlaps(plot_region) %>%
  trackCTSS(name = "CTSSs")

# Cluster track
cluster_track <- RSE %>%
  rowRanges %>%
  subsetByOverlaps(plot_region) %>%
  trackClusters(name = "Clusters", 
                col = NA, 
                showId = TRUE)

# Plot tracks together
plotTracks(list(axis_track, 
                ctss_track,
                cluster_track,
                tx_track),
           from = start(plot_region), 
           to = end(plot_region), 
           chromosome = as.character(seqnames(plot_region)))

#Location and expression of TSSs and enhancers
cluster_info <- RSE %>%
  rowData() %>%
  as.data.frame()

# Number of clusters
ggplot(cluster_info, aes(x = txType, fill = clusterType)) +
  geom_bar(alpha = 0.75, position = "dodge", color = "black") +
  scale_fill_colorblind("Cluster type") +
  labs(x = "Cluster annotation", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Expression of clusters
ggplot(cluster_info, aes(x = txType, 
                         y = log2(score / ncol(RSE)), 
                         fill = clusterType)) +
  geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.50, 0.75)) +
  scale_fill_colorblind("Cluster type") +
  labs(x = "Cluster annotation", y = "log2(TPM)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Analysing TSS shapes and sequences

# Select highly expressed TSSs
highTSSs <- subset(RSE, clusterType == 'TSS' & score / ncol(RSE) >= 10)

# Extended Data Figure 1d
# Calculate IQR as 10%-90% interval 
highTSSs <- calcShape(highTSSs, 
                      pooled = CTSSs, 
                      shapeFunction = shapeIQR, 
                      lower = 0.10, 
                      upper = 0.90)

highTSSs %>%
  rowData() %>%
  as.data.frame() %>%
  ggplot(aes(x = IQR)) +
  geom_histogram(binwidth = 1, 
                 fill = "#009E73", 
                 alpha = 0.75) +
  geom_vline(xintercept = 8, 
             linetype = "dashed", 
             alpha = 0.75, 
             color = "black") +
  facet_zoom(xlim = c(0,100)) +
  labs(x = "10-90% IQR", 
       y = "Frequency")

# Divide into groups
rowData(highTSSs)$shape <- ifelse(rowData(highTSSs)$IQR < 8, "Sharp", "Broad")

# Count group sizes
table(rowData(highTSSs)$shape)

# the core promoter sequences of the two classes of TSS candidates
promoter_seqs <- highTSSs %>%
  rowRanges() %>%
  swapRanges() %>%
  promoters(upstream = 40, downstream = 10) %>%
  getSeq(bsg, .)

promoter_seqs %>%
  as.character %>%
  split(rowData(highTSSs)$shape) %>%
  ggseqlogo(data = ., ncol = 2, nrow = 1) +
  theme_logo() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Finding candidates for interacting TSSs and enhancers

rowData(RSE)$clusterType <- RSE %>%
  rowData() %>%
  use_series("clusterType") %>%
  as_factor() %>%
  fct_relevel("TSS")

# Find links and calculate correlations
all_links <- RSE %>%
  swapRanges() %>%
  findLinks(maxDist = 1e5L,
            directional = "clusterType",
            inputAssay = "TPM",
            method = "kendall")

all_links


# Subset to only positive correlation
cor_links <- subset(all_links, estimate > 0 & p.value <0.05)

# Sort based on correlation
cor_links <- cor_links[order(cor_links$estimate, decreasing = TRUE)]

# Extract region around the link of interest
plot_region <- cor_links[1] %>% 
  boundingBox() %>% 
  linearize(1:2) %>%
  add(1000L)

# Cluster track
cluster_track <- RSE %>%
  subsetByOverlaps(plot_region) %>%
  trackClusters(name = "Clusters", 
                col = NA, 
                showId = TRUE)

# Link track
link_track <- cor_links %>%
  subsetByOverlaps(plot_region) %>%
  trackLinks(name = "Links",
             interaction.measure = "p.value",
             interaction.dimension.transform = "log",
             col.outside = "grey",
             plot.anchors = FALSE,
             col.interactions = "black")

# Plot tracks together
plotTracks(list(axis_track, 
                link_track,
                cluster_track,
                tx_track),
           from = start(plot_region), 
           to = end(plot_region), 
           chromosome = as.character(seqnames(plot_region)))

# Subset to only enhancers
Enhancers <- subset(RSE, clusterType == "Enhancer")

# Find stretches within 12.5 kbp
stretches <- findStretches(Enhancers, 
                           inputAssay = "TPM",
                           mergeDist = 12500L,
                           minSize = 5L,
                           method = "kendall")

stretches <- assignTxType(stretches, txModels = txdb)

# Sort by correlation
stretches <- stretches[order(stretches$aveCor, decreasing = TRUE)]
cor_links <- subset(stretches, aveCor > 0)

# Show the results
plot_region <- cor_links[1] + 1000

# Cluster track
cluster_track <- RSE %>%
  subsetByOverlaps(plot_region) %>%
  trackClusters(name = "Clusters", 
                col = NA, 
                showId = TRUE)

# CTSS track
ctss_track <- CTSSs %>%
  subsetByOverlaps(plot_region) %>%
  trackCTSS(name = "CTSSs")

# Stretch enhancer track
stretch_track <- stretches %>%
  subsetByOverlaps(plot_region) %>%
  AnnotationTrack(name = "Stretches", fill = "#009E73", col = NULL)

# Plot tracks together
plotTracks(list(axis_track, 
                stretch_track,
                cluster_track,
                ctss_track),
           from = start(plot_region), 
           to = end(plot_region), 
           chromosome = as.character(seqnames(plot_region)))

# Differential Expression analysis of TSSs, enhancers and genes

# Create DESeq2 object with blank design
dds_blind <- DESeqDataSet(RSE, design = ~ 1)

# Normalize and log transform
vst_blind <- vst(dds_blind, blind = TRUE)

##Figure 1a ###
plotPCA(vst_blind, "Class")

#p1 <- plotPCA(vst_blind, "Class")
#ggsave("plot/PCA for TSS.pdf", plot = p1, device = cairo_pdf, units = "in", dpi = 300)

# Specify design
dds <- DESeqDataSet(RSE, design = ~ Class)

# Fit DESeq2 model
dds <- DESeq(dds)

# Extract results
res <- results(dds,
               contrast = c("Class", "p16", "p12"),
               alpha = 0.05, 
               independentFiltering = TRUE, 
               tidy = TRUE) %>%
  bind_cols(as.data.frame(rowData(RSE))) %>%
  as_tibble()

res %>%
  top_n(n = -10, wt = padj) %>%
  dplyr::select(cluster = row, 
                clusterType, 
                txType, 
                baseMean, 
                log2FoldChange, 
                padj) %>%
  knitr::kable(caption = "Top differentially expressed TSS and enhancer candidates") 

ggplot(res, aes(x = log2(baseMean), 
                y = log2FoldChange, 
                color = padj < 0.05)) +
  geom_point(alpha = 0.25) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             alpha = 0.75) +
  facet_grid(clusterType ~ .)

#-------Extended Data Figure 1c--------------------------
sig_counts <- res %>%
  filter(padj < 0.05) %>%
  group_by(clusterType) %>%
  summarize(
    total_sig = n(),
    up = sum(log2FoldChange > 0, na.rm = TRUE),
    down = sum(log2FoldChange < 0, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(label = sprintf("Sig: %d\nUp: %d\nDown: %d", total_sig, up, down))

p <-ggplot(res, aes(x = log2(baseMean), 
                    y = log2FoldChange, 
                    color = padj < 0.05)) +
  geom_point(alpha = 0.25) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             alpha = 0.75) +
  facet_grid(clusterType ~ .) +
  geom_text(data = sig_counts,
            aes(x = Inf, y = -Inf, label = label), 
            inherit.aes = FALSE,                  
            color = "black",                      
            hjust = 1.1,                          
            vjust = -0.1,                          
            size = 4,                             
            fontface = "bold")

p_custom <- p + 
  theme(
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "#F0F0F0", color = NA),
    strip.text = element_text(color = "black") 
  )

print(p_custom)

#-----------------------------
#Extended Data Figure 4a

table(clusterType = rowRanges(RSE)$clusterType, 
      DE = res$padj < 0.05)

res_TSS <- res %>%
  dplyr::filter(clusterType == "TSS")

ggplot(res_TSS, aes(x = log2(baseMean), 
                    y = log2FoldChange, 
                    color = padj < 0.05)) +
  geom_point(alpha = 0.25) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             alpha = 0.75) +
  facet_grid(clusterType ~ .)

# Find top 10 DE enhancers
top10 <- res %>%
  filter(clusterType == "Enhancer", padj < 0.05) %>%
  group_by(log2FoldChange >= 0) %>%
  top_n(n = 5, wt = abs(log2FoldChange)) %>%
  pull(row)

# Extract expression values in tidy format
tidyEnhancers <- assay(RSE, "counts")[top10, ] %>%
  t() %>%
  as_tibble(rownames = "Sample") %>%
  mutate(class = RSE$Class) %>%
  gather(key = "Enhancer", 
         value = "Expression", 
         -Sample, -class, 
         factor_key = TRUE)

ggplot(tidyEnhancers, aes(x = class, 
                          y = Expression, 
                          fill = class)) +
  geom_dotplot(stackdir = "center", 
               binaxis = "y", 
               dotsize = 3) +
  facet_wrap(~ Enhancer, 
             ncol = 2, 
             scales = "free_y")

# Enrichment of DNA-binding motifs
cluster_seqs <- RSE %>% 
  rowRanges() %>%
  swapRanges() %>%
  unstrand() %>%
  add(500) %>% 
  getSeq(bsg, names = .)

# Extract motifs as PFMs
motif_pfms <- getMatrixSet(JASPAR2020, opts = list(species = "9606"))

# Look at the IDs and names of the first few motifs:
head(name(motif_pfms))

# Find matches
motif_hits <- matchMotifs(motif_pfms, subject = cluster_seqs)

# Matches are returned as a sparse matrix:
motifMatches(motif_hits)[1:5, 1:5]

table(TF = motifMatches(motif_hits)[,"MA0636.1"],
      DE = res$padj < 0.05) %>%
  print() %>%
  fisher.test()

# ---------- Gene-level differential expression ---------------
RSE <- assignGeneID(RSE, geneModels = txdb)

GSE <- RSE %>%
  subset(clusterType == "TSS") %>%
  quantifyGenes(genes = "geneID", inputAssay = "counts")

rowRanges(GSE["259307",])
print(rownames(GSE)[2])
# Translate symbols
rowData(GSE)$symbol <- mapIds(odb, 
                              keys = rownames(GSE), 
                              column = "SYMBOL", 
                              keytype = "ENTREZID")

# Create DGElist object
dge <- DGEList(counts = assay(GSE, "counts"),
               genes = as.data.frame(rowData(GSE)))

# Calculate normalization factors
dge <- calcNormFactors(dge)

mod <- model.matrix(~ Class - 1, data = colData(GSE))

contr.matrix <- makeContrasts(p16vsp3 = Classp16-Classp03,
                              p14vsp3 = Classp14-Classp03, 
                              p9vsp3 = Classp09-Classp03,
                              p14vsp9 = Classp14-Classp09,
                              p12vsp9 = Classp12-Classp09,
                              p14vsp12 = Classp14-Classp12,
                              p16vsp14 = Classp16-Classp14,
                              p18vsp16 = Classp18-Classp16,
                              p16vsp12 = Classp16-Classp12,
                              p18vsp3 = Classp18-Classp03,
                              levels = colnames(mod))


# Model mean-variance using voom
v <- voom(dge, design = mod,plot=TRUE)

# Fit and shrink DE model
vfit <- lmFit(v, design = mod)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
eb <- eBayes(vfit, robust = TRUE)
express <- eb$coefficients[1]
print(colnames(eb$coefficients))
plotSA(eb, main="Final model: Mean-variance trend")

# Summarize the results
dt <- decideTests(eb)

# Global summary
dt %>% 
  summary() %>% 
  knitr::kable(caption = "Global summary of differentially expressed genes.")

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))

result_table <-  topTable(eb, coef = "p18vsp3", n = Inf) %>%
  dplyr::select(symbol, nClusters, AveExpr, logFC, adj.P.Val) %>%
  knitr::kable(caption = "Top differentially expressed genes.")
head(result_table)

#write.csv(result_table, "C:/Users/Jiaqi/Documents/Code/cage/Output/p16vsp3.csv", row.names = FALSE)


GO <- goana(eb, coef = "p18vsp3", species = "Hs", trend = TRUE)

topGO(GO, ontology = "BP", number = 20) %>%
  knitr::kable(caption = "Top enriched or depleted GO-terms.")

KEGG <- kegga(eb, coef = "p16vsp3", species = "Hs", trend = TRUE)

# Show top hits
topKEGG(KEGG, number = 20) %>%
  knitr::kable(caption = "Top enriched or depleted KEGG-terms.")

#  PLOT DG

p16vsp12_results <- topTable(eb, coef = "p18vsp3", number = Inf, sort.by = "P")
# Screen DE
sig_genes <- p16vsp12_results %>%
  filter(adj.P.Val < 0.05 & logFC > 0.6) %>%  
  pull(symbol)

sig_genes_entrez <- bitr(sig_genes, fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = odb)

# GO enrichment
ego <- enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = odb, 
                keyType = "ENTREZID", ont = "BP", 
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(ego, showCategory = 10) + ggtitle("GO Enrichment for p16 vs p12")

# KEGG enrichment
ekegg <- enrichKEGG(gene = sig_genes_entrez$ENTREZID, organism = "hsa", 
                    pAdjustMethod = "BH", pvalueCutoff = 0.05)
dotplot(ekegg, showCategory = 10) + ggtitle("KEGG Enrichment for p16 vs p12")

#  results 
volcano_data <- p16vsp12_results %>%
  mutate(
    diffexpressed = case_when(
      logFC > 0.6 & adj.P.Val < 0.05 ~ "Upregulated",
      logFC < -0.6 & adj.P.Val < 0.05 ~ "Downregulated",
      TRUE ~ "NO"
    ),
    importance = abs(logFC) * -log10(adj.P.Val) 
  )

# DE
filtered_data <- volcano_data %>%
  filter(diffexpressed != "NO")
filtered_data
# top15 DE
top_genes <- filtered_data %>%
  arrange(desc(importance)) %>%
  head(15) %>%
  pull(symbol)

# Add TAG
volcano_data$delabel <- NA
volcano_data$delabel[volcano_data$symbol %in% top_genes] <- volcano_data$symbol[volcano_data$symbol %in% top_genes]

summary(volcano_data$gene)
summary(volcano_data$symbol)


# Biostatsquid theme
theme_set(theme_classic(base_size = 12) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

ggplot(data=volcano_data, aes(x = logFC, y = -log10(adj.P.Val), col = diffexpressed, label = delabel))+
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed')+
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed')+
  geom_point(size = 2)+
  scale_color_manual(values = c("#00AFBB", "grey", "#F8766D"), 
                     labels = c("Downregulated", "Not significant", "Upregulated"))+
  labs(color = 'Differential expressed genes', 
       x = expression("logFC"), y = expression("-log"[10]*"adj.P.Val"))+
  scale_x_continuous(breaks = seq(-10, 10, 2)) + 
  ggtitle('HUVEC p16 vs p12') + 
  geom_text_repel(max.overlaps = Inf) 


##Figure 1b,c###
##--------------------------GSEA FOR SEN-MAYO------------------------------##
#prepare the data for analysis
library(fgsea)
library(msigdbr)
library(GSEABase)

group_info <- data.frame(
  p03 = c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  p09 = c(0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  p12 = c(0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  p14 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0),
  p16 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0),
  p18 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1)
)

eset<-dge
eset$samples$group <- apply(group_info, 1, function(x) {
  group_names <- colnames(group_info)
  paste(group_names[x == 1], collapse = ",")
})
eset$counts <- normalizeBetweenArrays(log2(eset$counts+1), method="quantile")
eset$counts <- eset$counts[order(rowMeans(eset$counts), decreasing = TRUE), ]
eset$counts <- eset$counts[1:12000, ]
head(eset$counts)
pathwaysDF <- msigdbr(species = "Homo sapiens",  category="H")
head(pathwaysDF)
pathways <- split(as.character(pathwaysDF$entrez_gene), pathwaysDF$gs_name)
head(pathways)

#matching_pathways <- grep("SAUL_SEN_MAYO", names(pathways), value = TRUE)
#print(matching_pathways)

gmt.file <- "C:/Users/Jiaqi/Documents/Code/R/aging/SAUL_SEN_MAYO.v2023.2.Hs.gmt"
sen <- gmtPathways(gmt.file)
sen <- bitr(sen$SAUL_SEN_MAYO, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

sen <- as.character(sen$ENTREZID)
sen

set.seed(1)
gesecaRes <- geseca(pathways, eset$counts, minSize = 15, maxSize = 500, eps = 0)
head(gesecaRes, 10)

plotCoregulationProfile(pathway=sen, 
                        E=eset$counts, conditions=eset$samples$group)
plotGesecaTable(gesecaRes |> head(10), pathways, E=eset$counts,
                colors = c("#00AFBB", "#FFFFFF", "#F8766D"))

##Figure 1d###
##---------------------------------mfuzz-------------------------------------##
mfuzz_matrix <- eset$counts
group_info <- eset$samples$group
mset <- t(apply(mfuzz_matrix, 1, function(x) tapply(x, group_info, mean)))
head(mset)
mset <- as.matrix(mset)
gene_covert<-rownames(mset)

gene_covert <- bitr(gene_covert, fromType = "ENTREZID", 
                    toType = "SYMBOL", 
                    OrgDb = org.Hs.eg.db)

gene_ids <- rownames(mset)

gene_symbols <- gene_covert$SYMBOL[match(gene_ids, gene_covert$ENTREZID)]

# create matrix with gene symbol
mset <- data.frame(gene = gene_symbols, mset,stringsAsFactors = FALSE)


names(mset)[1] <- ""
any(is.na(mset[, 1]))
mset[is.na(mset[, 1]), 1] <- "NA"
rownames(mset) <- mset[, 1]
mset <- mset[, -1]

mset <- new("ExpressionSet", exprs = as.matrix(mset))
mset<-filter.NA(mset,thres = 0.25)
mset<- fill.NA(mset,mode="Wknn")
mset<-filter.std(mset, min.std = 0.1)
mset <- standardise(mset)
set.seed(123)
m1 <- mestimate(mset)
m1
cl <- mfuzz(mset,c=6,m=m1)
mfuzz.plot2(mset,cl=cl,mfrow=c(2,3),time.labels=c("p03","p09","p12","p14","p16","p18"))

tmp<-acore(mset,cl,min.acore = 0.6)
tail(tmp)
#up_gene <- tmp[[3]]
#write.csv(up_gene, "C:/Users/Jiaqi/Documents/Code/R/aging/output/upgene.csv", row.names = FALSE)

#----------------------------------------------------------------
# Enrichment
# 1. Extract SYMBOL of every clusters and name as list 

gene_list <- lapply(tmp, function(cluster_df) {
  cluster_df$NAME  # extrac gene name
})

names(gene_list) <- paste0(seq_along(gene_list))

# check results
str(gene_list)

# convert SYMBOL to ENTREZID，keeo cluster names
entrez_list <- lapply(gene_list, function(symbols) {
  mapped <- bitr(
    geneID = symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = odb
  )
  # extract ENTREZID and convert to character
  if (nrow(mapped) == 0) {
    return(character(0))  
  }
  as.character(unique(mapped$ENTREZID))
})

str(entrez_list)

ego <- compareCluster(
  geneCluster = entrez_list,
  fun = "enrichGO",
  OrgDb = odb,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# plot
dotplot(ego, showCategory = 3) +
  labs(title = "GO Enrichment")



cnetplot(ego)

ck <- compareCluster(
  geneCluster = entrez_list,
  fun = "enrichKEGG",
  organism="hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

dotplot(ck, showCategory = 4) +
  labs(title = "KEGG Enrichment")

cnetplot(ck,
         showCategory = 5,
         layout = "kk",
         circular = FALSE,
         colorEdge = TRUE,
         cex_label_category = 0.8,
         cex_label_gene = 0.6)

##----------------------------------------------------------------##
##----------------------------------------------------------------##

#Figure 2a-d# Extended Data Figure 1e, f

# packages
library(DiffLogo)
library(Biostrings)
library(ggseqlogo)
library(ggplot2)
library(tidyverse)
library(Biostrings)

#---------------------------------

# Extract results
res <- results(dds,
               contrast = c("Class", "p09", "p03"),
               alpha = 0.05, 
               independentFiltering = TRUE, 
               tidy = TRUE) %>%
  bind_cols(as.data.frame(rowData(RSE))) %>%
  as_tibble()

res %>%
  top_n(n = -10, wt = padj) %>%
  dplyr::select(cluster = row, 
                clusterType, 
                txType, 
                baseMean, 
                log2FoldChange, 
                padj) %>%
  knitr::kable(caption = "Top differentially expressed TSS and enhancer candidates") 

tss_ovy <- res

res_ovy <- subset(tss_ovy, clusterType == "TSS")
res_up <- subset(res_ovy, padj < 0.05 & log2FoldChange > 0) 
res_down<-subset(res_ovy, padj < 0.05 & log2FoldChange < 0) 

# Select highly expressed TSSs
highTSSs <- subset(RSE, clusterType == 'TSS' & score / ncol(RSE) >= 10)

# Calculate IQR as 10%-90% interval 
highTSSs <- calcShape(highTSSs, 
                      pooled = CTSSs, 
                      shapeFunction = shapeIQR, 
                      lower = 0.10, 
                      upper = 0.90)

highTSSs %>%
  rowData() %>%
  as.data.frame() %>%
  ggplot(aes(x = IQR)) +
  geom_histogram(binwidth = 1, 
                 fill = "hotpink", 
                 alpha = 0.75) +
  geom_vline(xintercept = 8, 
             linetype = "dashed", 
             alpha = 0.75, 
             color = "black") +
  facet_zoom(xlim = c(0,100)) +
  labs(x = "10-90% IQR", 
       y = "Frequency")

# Divide into groups
rowData(highTSSs)$shape <- ifelse(rowData(highTSSs)$IQR < 8, "Sharp", "Broad")

# Count group sizes
table(rowData(highTSSs)$shape)

avg_width <- highTSSs %>%
  rowData() %>%
  as.data.frame() %>%
  group_by(shape) %>%
  summarise(avg_width = mean(IQR))
avg_width

y_index <- res_down$row
s_index <- res_up$row
s_highTSSs <- highTSSs[rownames(highTSSs) %in% s_index, ]
y_highTSSs<-highTSSs[rownames(highTSSs) %in% y_index , ]

sen_seqs <- s_highTSSs %>%
  rowRanges() %>%
  swapRanges() %>%
  promoters(upstream = 40, downstream = 10) %>%
  getSeq(bsg, .)

sen_seqs %>%
  as.character %>%
  split(rowData(s_highTSSs)$shape) %>%
  ggseqlogo(data = ., ncol = 2, nrow = 1) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme_logo() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

young_seqs <- y_highTSSs %>%
  rowRanges() %>%
  swapRanges() %>%
  promoters(upstream = 40, downstream = 10) %>%
  getSeq(bsg, .)

young_seqs %>%
  as.character %>%
  split(rowData(y_highTSSs)$shape) %>%
  ggseqlogo(data = ., ncol = 2, nrow = 1) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme_logo() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#--------------------------------

# Define a function
get_global_freq <- function(seqs, group_name, passage_name) {
  # calculate base 
  alphabet_counts <- alphabetFrequency(seqs, baseOnly = TRUE)
  total_counts <- colSums(alphabet_counts[, c("A", "C", "G", "T")])
  
  # convert to frequency
  freqs <- total_counts / sum(total_counts)
  
  # concert to data.frame 
  data.frame(
    Base = names(freqs),
    Frequency = as.numeric(freqs),
    Shape = group_name,
    Passage = passage_name
  )
}

# ------

# 1. extract Sharp group
sharp_up_freq <- get_global_freq(sen_seqs[rowData(s_highTSSs)$shape == "Sharp"], "Sharp", "p09_UP")
sharp_down_freq <- get_global_freq(young_seqs[rowData(y_highTSSs)$shape == "Sharp"], "Sharp", "p03_UP")

# 2. extract Sharp group
broad_up_freq <- get_global_freq(sen_seqs[rowData(s_highTSSs)$shape == "Broad"], "Broad", "p09_UP")
broad_down_freq <- get_global_freq(young_seqs[rowData(y_highTSSs)$shape == "Broad"], "Broad", "p03_UP")

# 3. combine 
plot_data <- bind_rows(sharp_up_freq, sharp_down_freq, broad_up_freq, broad_down_freq)


ggplot(plot_data, aes(x = Base, y = Frequency, fill = Passage)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single"), color = "black", width = 0.7) +
  facet_wrap(~Shape) + 
  scale_fill_brewer(palette = "Set1") + 
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 0.5)) +
  labs(
    title = "Base Composition Preference across Passages",
    subtitle = "Comparing global nucleotide usage in Sharp vs Broad TSS",
    x = "Nucleotide",
    y = "Average Frequency (%)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "top"
  )

plot_data %>%
  group_by(Shape, Passage) %>%
  summarise(GC_Content = sum(Frequency[Base %in% c("G", "C")]))



# ------

up_sharp_seqs_fix   <- sen_seqs[rowData(s_highTSSs)$shape == "Sharp"]     
down_sharp_seqs_fix <- young_seqs[rowData(y_highTSSs)$shape == "Sharp"]   

get_freq <- function(x) {
  counts <- alphabetFrequency(x, baseOnly = TRUE)
  total <- colSums(counts[, c("A", "C", "G", "T")])
  total / sum(total)
}

freq_p09 <- get_freq(up_sharp_seqs_fix) 
freq_p03 <- get_freq(down_sharp_seqs_fix)

plot_data_sharp <- data.frame(
  Base = c("A", "C", "G", "T"),
  Ratio = as.numeric(freq_p09 / freq_p03),
  Passage = "P9"
)

baseline_data <- data.frame(
  Base = c("A", "C", "G", "T"),
  Ratio = 1,
  Passage = "P3"
)

final_df <- bind_rows(baseline_data, plot_data_sharp)

ggplot(final_df, aes(x = Base, y = Ratio, fill = Passage)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.7), 
           color = "black", 
           width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgrey", size = 0.8) +
  scale_fill_manual(values = c("P3" = "#00AFBB", 
                               "P9" = "#F8766D")) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, 1.3), 
                     breaks = seq(0, 1.3, 0.1)) +
  labs(
    title = "Relative change in base frequency (P9/P3 ratio)",
    x = "Nucleotide",
    y = "Relative Usage"
  ) +
  theme_classic() + 
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )
#---------------------------
# Extract results
res <- results(dds,
               contrast = c("Class", "p16", "p12"),
               alpha = 0.05, 
               independentFiltering = TRUE, 
               tidy = TRUE) %>%
  bind_cols(as.data.frame(rowData(RSE))) %>%
  as_tibble()

res %>%
  top_n(n = -10, wt = padj) %>%
  dplyr::select(cluster = row, 
                clusterType, 
                txType, 
                baseMean, 
                log2FoldChange, 
                padj) %>%
  knitr::kable(caption = "Top differentially expressed TSS and enhancer candidates") 

tss_ovy <- res

res_ovy <- subset(tss_ovy, clusterType == "TSS")
res_up <- subset(res_ovy, padj < 0.05 & log2FoldChange > 0) 
res_down<-subset(res_ovy, padj < 0.05 & log2FoldChange < 0) 

# Select highly expressed TSSs
highTSSs <- subset(RSE, clusterType == 'TSS' & score / ncol(RSE) >= 10)

# Calculate IQR as 10%-90% interval 
highTSSs <- calcShape(highTSSs, 
                      pooled = CTSSs, 
                      shapeFunction = shapeIQR, 
                      lower = 0.10, 
                      upper = 0.90)

highTSSs %>%
  rowData() %>%
  as.data.frame() %>%
  ggplot(aes(x = IQR)) +
  geom_histogram(binwidth = 1, 
                 fill = "hotpink", 
                 alpha = 0.75) +
  geom_vline(xintercept = 8, 
             linetype = "dashed", 
             alpha = 0.75, 
             color = "black") +
  facet_zoom(xlim = c(0,100)) +
  labs(x = "10-90% IQR", 
       y = "Frequency")

# Divide into groups
rowData(highTSSs)$shape <- ifelse(rowData(highTSSs)$IQR < 8, "Sharp", "Broad")

# Count group sizes
table(rowData(highTSSs)$shape)

avg_width <- highTSSs %>%
  rowData() %>%
  as.data.frame() %>%
  group_by(shape) %>%
  summarise(avg_width = mean(IQR))
avg_width

y_index <- res_down$row
s_index <- res_up$row
s_highTSSs <- highTSSs[rownames(highTSSs) %in% s_index, ]
y_highTSSs<-highTSSs[rownames(highTSSs) %in% y_index , ]

sen_seqs <- s_highTSSs %>%
  rowRanges() %>%
  swapRanges() %>%
  promoters(upstream = 40, downstream = 10) %>%
  getSeq(bsg, .)

sen_seqs %>%
  as.character %>%
  split(rowData(s_highTSSs)$shape) %>%
  ggseqlogo(data = ., ncol = 2, nrow = 1) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme_logo() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

young_seqs <- y_highTSSs %>%
  rowRanges() %>%
  swapRanges() %>%
  promoters(upstream = 40, downstream = 10) %>%
  getSeq(bsg, .)

young_seqs %>%
  as.character %>%
  split(rowData(y_highTSSs)$shape) %>%
  ggseqlogo(data = ., ncol = 2, nrow = 1) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme_logo() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#--------------------------------
sharp_up_freq <- get_global_freq(sen_seqs[rowData(s_highTSSs)$shape == "Sharp"], "Sharp", "p09_UP")
sharp_down_freq <- get_global_freq(young_seqs[rowData(y_highTSSs)$shape == "Sharp"], "Sharp", "p03_UP")
broad_up_freq <- get_global_freq(sen_seqs[rowData(s_highTSSs)$shape == "Broad"], "Broad", "p09_UP")
broad_down_freq <- get_global_freq(young_seqs[rowData(y_highTSSs)$shape == "Broad"], "Broad", "p03_UP")

plot_data <- bind_rows(sharp_up_freq, sharp_down_freq, broad_up_freq, broad_down_freq)


ggplot(plot_data, aes(x = Base, y = Frequency, fill = Passage)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single"), color = "black", width = 0.7) +
  facet_wrap(~Shape) + 
  scale_fill_brewer(palette = "Set1") + 
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 0.5)) +
  labs(
    title = "Base Composition Preference across Passages",
    subtitle = "Comparing global nucleotide usage in Sharp vs Broad TSS",
    x = "Nucleotide",
    y = "Average Frequency (%)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "top"
  )

plot_data %>%
  group_by(Shape, Passage) %>%
  summarise(GC_Content = sum(Frequency[Base %in% c("G", "C")]))


up_sharp_seqs_fix   <- sen_seqs[rowData(s_highTSSs)$shape == "Sharp"]     
down_sharp_seqs_fix <- young_seqs[rowData(y_highTSSs)$shape == "Sharp"]   

get_freq <- function(x) {
  counts <- alphabetFrequency(x, baseOnly = TRUE)
  total <- colSums(counts[, c("A", "C", "G", "T")])
  total / sum(total)
}

freq_p16 <- get_freq(up_sharp_seqs_fix) 
freq_p12 <- get_freq(down_sharp_seqs_fix)

plot_data_sharp <- data.frame(
  Base = c("A", "C", "G", "T"),
  Ratio = as.numeric(freq_p16 / freq_p12),
  Passage = "P16"
)

baseline_data <- data.frame(
  Base = c("A", "C", "G", "T"),
  Ratio = 1,
  Passage = "P12"
)

final_df <- bind_rows(baseline_data, plot_data_sharp)

ggplot(final_df, aes(x = Base, y = Ratio, fill = Passage)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.7), 
           color = "black", 
           width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgrey", size = 0.8) +
  scale_fill_manual(values = c("P12" = "#00AFBB", 
                               "P16" = "#F8766D")) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, 1.3), 
                     breaks = seq(0, 1.3, 0.1)) +
  labs(
    title = "Relative change in base frequency (P16/P12 ratio)",
    x = "Nucleotide",
    y = "Relative Usage"
  ) +
  theme_classic() + 
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )


## Figure 6d ##
#----------------------------------------------------------
library(GenomicRanges)
library(Gviz)
library(CAGEfightR)

# 1. target range
plot_region <- GRanges("chr12", IRanges(53380000, 53380950))

# 2. extract TSS 
CTSSs_sub <- subsetByOverlaps(CTSSs, plot_region)

# 3. Separate group
groups <- c("p03", "p09", "p12", "p16")
CTSSs_sub <- CTSSs_sub[, CTSSs_sub$Class %in% groups]

max_tpm <- max(assay(CTSSs_sub, "TPM"), na.rm = TRUE)
if(is.infinite(max_tpm) || max_tpm == 0) max_tpm <- 1 
ylim_range <- c(-max_tpm * 1.05, max_tpm * 1.05)
tracks_list <- list()

for (grp in groups) {
  
  ctss_grp <- CTSSs_sub[, CTSSs_sub$Class == grp]
  ctss_grp <- calcPooled(ctss_grp) 
  
  track_name <- paste0(grp, "\n")
  
  dt <- trackCTSS(ctss_grp, name = track_name)
  
  displayPars(dt) <- list(
    ylim = ylim_range,
    rotation.title = 0,      
    just.title = "right",                 
    background.title = "transparent",     
    col.border.title = "transparent",     
    col.title = "black",                   
    cex.title = 1.0,         
    cex.axis = 0.8,          
    lwd = 1.5,
    
    col = c("#00BFC4", "#F8766D")
  )
  
  tracks_list[[grp]] <- dt
}


# 
axis_track <- GenomeAxisTrack(fontsize = 12) 

displayPars(axis_track) <- list(
  labelPos = "alternating",   
  cex.axis = 0.8,             
  fontfamily = "Arial",       
  fontcolor = "grey",         
  col.axis = "grey",          
  distFromAxis = 1            
)

#  Arial
displayPars(tx_track) <- list(
  rotation.title = 0, 
  just.title = "right",
  background.title = "transparent",
  col.border.title = "transparent",
  col.title = "black",         
  cex.title = 1.0, 
  fontfamily.title = "Arial", 
  fontfamily.group = "Arial",
  fontsize.group = 12 
)

#  ordering
ordered_tracks <- tracks_list[order(names(tracks_list))]

#  Plotting
plotTracks(
  c(list(axis_track), ordered_tracks, list(tx_track)),
  from = start(plot_region),
  to = end(plot_region),
  chromosome = as.character(seqnames(plot_region)),
  title.width = 0.6,      
  margin = 10,             
  innerMargin = 4,         
  fontsize = 16,           
  sizes = c(1.0, rep(1, length(ordered_tracks)), 1.0) 
)
