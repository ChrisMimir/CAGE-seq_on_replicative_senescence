 # CAGEfightR data to ANANSE #

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")

BiocManager::install("CAGEfightR")
BiocManager::install(c("rtracklayer","InteractionSet","stringr",
                       "dplyr","tidyr","DESeq2"))

BiocManager::install(c("TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "org.Hs.eg.db"))
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")


#remove.packages(row.names(installed.packages()))

knitr::opts_chunk$set(echo = TRUE)

library(CAGEfightR)
library(rtracklayer)
library(InteractionSet)
library(stringr)
library(dplyr)
library(tidyr)
library(DESeq2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  
library(BSgenome.Hsapiens.UCSC.hg38)  
library(org.Hs.eg.db)  

# Genome info 
UCSCGenome <- "hg38" 
BS_genome <- "BSgenome.Hsapiens.UCSC.hg38"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 
odb <- org.Hs.eg.db


# Data
setwd("C:/Users/Jiaqi/Documents/Code/R/ANANSE-CAGE/work")
work_dir <- getwd()

data = "C:/Users/Jiaqi/Documents/Code/R/ANANSE-CAGE/data"
print(list.files(path = data))

# To convert ctss files to bigwig
PREPROCESSING = FALSE

# Select samples -  Change accordingly!
source_type <- "p9"
target_types <- "p16"

UNEXPRESSED <- 0  # minimal TPM

if (PREPROCESSING == TRUE){
  hg38 = SeqinfoForUCSCGenome("hg38")
  bedfiles <- list.files(path = data, pattern = "*.bed", full.names = TRUE)
  if (length(bedfiles) == 0){print("No BED files in data directory. Pre-processing failed.")}
  
  for (i in bedfiles) {
    print(i)
    
    name_plus <- str_replace(i, ".bed", ".CTSS.raw.plus.bw")
    name_minus <- str_replace(i, ".bed", ".CTSS.raw.minus.bw")
    convertBED2BigWig(input=i,
                      outputPlus=name_plus,
                      outputMinus=name_minus,
                      genome=hg38)
  }
} else {
  print("Pre-processing skipped.")
}

for (type in target_types){
  out_dir <- paste(source_type, "_to_", type, sep = "")
  dir.create(file.path(work_dir, out_dir))} 
target_type = type

genomeInfo <- seqinfo(txdb)

# Create two named BigWigFileList-objects 
pattern_plus <- paste("(", source_type, "|", target_type, ").*plus",
                      sep = "")
pattern_minus <- paste("(", source_type, "|", target_type, ").*minus",
                       sep = "")

bw_plus <- list.files(path = data, pattern = pattern_plus, full.names = T)
bw_minus <- list.files(path = data, pattern = pattern_minus, full.names = T)

# reorder files to source_type first and target_type second
if (grepl(source_type, basename(bw_plus[1])) == F){
  print("WARNING: Reorder file order.")
  bw_plus <- sort(bw_plus, decreasing = T)}
if (grepl(source_type, basename(bw_minus[1])) == F){
  bw_minus <- sort(bw_minus, decreasing = T)} 

bw_plus <- BigWigFileList(bw_plus)
bw_minus <- BigWigFileList(bw_minus)

names(bw_plus) <- sub(".CTSS.raw.plus.bw", "", basename(bw_plus))
names(bw_minus) <- sub(".CTSS.raw.minus.bw", "", basename(bw_minus))

source_number <- length(grep(source_type, names(bw_plus), value = F))
target_number <- length(grep(target_type, names(bw_plus), value = F))

# Set minimum samples manually
minimum_samples <- 2

# Sanity check
print(names(bw_plus))
print(names(bw_minus))

# Quantify CTSS, normalize, and pool CTSSs
CTSSs <- quantifyCTSSs(plusStrand = bw_plus, 
                       minusStrand = bw_minus, 
                       genome = genomeInfo)

CTSSs <- calcTPM(CTSSs, inputAssay="counts", outputAssay="TPM", 
                 outputColumn = "totalTags")

CTSSs <- calcPooled(CTSSs, inputAssay = "TPM")

# Remove excess noise
# (Default) Count number of samples with MORE ( > ) than 0 counts:
CTSSs <- calcSupport(CTSSs, 
                     inputAssay="counts", 
                     outputColumn="support", 
                     unexpressed=0)
table(rowRanges(CTSSs)$support)
# Adjust "support" to discard CTSSs only in 'n' samples
# Default = 1
# table(rowRanges(CTSSs)$support)
supportedCTSSs <- subset(CTSSs, support > 2)
supportedCTSSs <- calcTPM(supportedCTSSs, totalTags="totalTags")
supportedCTSSs <- calcPooled(supportedCTSSs)


# Unidirectional clustering
prefiltered_TCs <- clusterUnidirectionally(supportedCTSSs, 
                                           pooledCutoff=0, 
                                           mergeDist=20)

Unidirectional <- quantifyClusters(CTSSs,
                                   clusters=prefiltered_TCs,
                                   inputAssay="counts")

Unidirectional <- calcTPM(Unidirectional, 
                          totalTags = "totalTags")

# Only TSSs expressed at more than 1 TPM in more than 'n' samples
# Default: unexpressed = 0, minSamples = 1
Unidirectional <- subsetBySupport(Unidirectional,
                                  inputAssay="TPM",
                                  unexpressed=0,
                                  minSamples=1)

# Annotation (Transcript Models)
Unidirectional <- assignTxID(Unidirectional,
                             txModels=txdb,
                             outputColumn="txID")

Unidirectional <- assignTxType(Unidirectional,
                               txModels=txdb,
                               outputColumn="txTyp")

Unidirectional <- assignTxType(Unidirectional,
                               txModels=txdb,
                               outputColumn="peakTxType",
                               swap="thick")

Unidirectional <- assignGeneID(Unidirectional,
                               geneModels=txdb,
                               outputColumn="geneID")

# Match IDs to symbols
symbols_uni <- mapIds(odb,
                      keys=rowRanges(Unidirectional)$geneID,
                      keytype="ENTREZID",
                      column="SYMBOL")

# Add to object
rowRanges(Unidirectional)$symbol <- as.character(symbols_uni)
Unidirectional <- assignMissingID(Unidirectional,
                                  outputColumn="symbol")

Unidirectional <- assignGeneID(Unidirectional,
                               geneModels=txdb,
                               outputColumn="geneID")

# Bidirectional clustering
# Default: balanceThreshold=0.95
BCs <- clusterBidirectionally(CTSSs, balanceThreshold=0.95)
BCs <- calcBidirectionality(BCs, samples=CTSSs)

# Remove (excess noise) BCs that are not observed in > 'n' samples
# Default = 0
# table(BCs$bidirectionality)
BCs <- subset(BCs, bidirectionality > 0)

Bidirectional <- quantifyClusters(CTSSs,
                                  clusters=BCs,
                                  inputAssay="counts")

Bidirectional <- calcTPM(Bidirectional, 
                         totalTags = "totalTags")

# Remove excess noise: Only keep BCs expressed in more than 'n' samples
# Default minSamples = 1, unexpressed = 0
Bidirectional <- subsetBySupport(Bidirectional,
                                 inputAssay="counts",
                                 unexpressed=0,
                                 minSamples=1)


# Genelevel
genelevel <- quantifyGenes(Unidirectional,
                           genes="symbol",
                           inputAssay="counts")
genelevel$group <- factor(c(rep(source_type, source_number), 
                            rep(target_type, target_number)))
# genelevel$group <- relevel(genelevel$group, target_type)
print(genelevel$group)


# --------------------------------- DEseq2 --------------------------------- # 

dds_genelevel <- DESeqDataSet(genelevel, ~group)
dds <- DESeq(dds_genelevel)
res <- results(dds, contrast = c("group", target_type, source_type))


# --------------------------------- EXPORT --------------------------------- #      

# Enhancers & Regions
enhancers <- assays(Bidirectional)$TPM
print(colnames(enhancers))

# Enhancers source
average_tpm_source_enhancer <- rowMeans(enhancers[,c(1:source_number)])
source_enhancers <- data.frame(rownames(enhancers), 
                               Means = average_tpm_source_enhancer)

rownames(source_enhancers) <- NULL
colnames(source_enhancers) <- c("pos", "CAGE")
source_enhancers_name <- paste(out_dir, "enhancers_source.tsv", sep = "/")
write.table(source_enhancers, 
            source_enhancers_name,
            row.names = F, 
            quote = FALSE, 
            sep = "\t")

# Enhancers target 
average_tpm_target_enhancer <- rowMeans(enhancers[,c((source_number + 1):(source_number + target_number))])
target_enhancers <- data.frame(rownames(enhancers), 
                               Means = average_tpm_target_enhancer)
rownames(target_enhancers) <- NULL
colnames(target_enhancers) <- c("pos", "CAGE")
target_enhancers_names <- paste(out_dir, "enhancers_target.tsv", sep = "/")
write.table(target_enhancers, 
            target_enhancers_names, 
            row.names = F, 
            quote = FALSE, 
            sep = "\t")


# Gene expression (TPM)
genes_tpm <- calcTPM(genelevel)
genes_tpm <- assays(genes_tpm)$TPM
row.names.remove <- c("NULL")
genes_tpm <- genes_tpm[!(row.names(genes_tpm) %in% row.names.remove),]

for (i in colnames(genes_tpm)) {
  print(i)
  type_tpm <- data.frame(genes_tpm[,i])
  type_tpm <- tibble::rownames_to_column(type_tpm, "resid")
  colnames(type_tpm) <- c("resid", "tpm")
  
  if (grepl(source_type, i) == T){
    file_name <- str_replace(i, source_type, "")
    file_name <- paste(out_dir, "/TPM_source", file_name, ".txt", sep="")
    write.table(type_tpm, file_name, row.names = F, quote = FALSE, sep = "\t")
  } else {
    file_name <- str_replace(i, target_type, "")
    file_name <- paste(out_dir, "/TPM_target", file_name, ".txt", sep="")
    write.table(type_tpm, file_name, row.names = F, quote = FALSE, sep = "\t")
  }
}


# Differentially expressed genes
de_genes <- res[,c(2,6)]
row.names.remove <- c("NULL")
de_genes <- de_genes[!(row.names(de_genes) %in% row.names.remove),]
de_genes <- as.data.frame(de_genes)
de_genes <- tibble::rownames_to_column(de_genes, "resid")
de_genes <- drop_na(de_genes)

de_genes_name <- paste(out_dir, "DE_genes.txt", sep = "/")
write.table(de_genes, de_genes_name, row.names = F, quote = FALSE, sep = "\t") 
