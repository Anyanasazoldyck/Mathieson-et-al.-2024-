library(DESeq2)
library(pheatmap)


sc_data<- readRDS("data/integrated_cafs_seurat.rds")
ncol(sc_data) #110 cells

Idents(sc_data) <- sc_data$subtype


cts <- AggregateExpression(sc_data, 
                           group.by = c("subtype", "patient_id"),
                           assays = 'SCT',
                           return.seurat = FALSE)

cts <- cts$SCT

# convert to data.frame
cts <- as.data.frame(cts)
cts = subset(cts,apply(cts, 1, mean) >= 1)


#make a sample sheet
sample_info <- data.frame(

  sample = colnames(cts),
  
  condition = ifelse(grepl("CAF-S1", colnames(cts)), "CAF_S1", "CAF_S5")
)

sample_info$condition <- factor(sample_info$condition)

# RUNDESEQ2####
dds <- DESeqDataSetFromMatrix(
  countData = round(cts),  # DESeq2 needs integer counts
  colData = sample_info,
  design = ~ condition
)
dds = DESeq(dds)
em = as.data.frame(counts(dds, normalized=TRUE))
em = round(em,2)
em_out = cbind(row.names(em), em)
names(em_out)[1] = 'ID'
write.table(em_out, file="result/EM.csv",row.names=FALSE, sep="\t", quote = FALSE)


#DE ::::S5_vs_S1
S5_vs_S1 <- as.data.frame(results(dds,c("condition","CAF_S5","CAF_S1")))
S5_vs_S1$mlogp <- -log10(S5_vs_S1$padj)
S5_vs_S1$sig <- ifelse(S5_vs_S1$padj <0.05,"significant","nonsignificant")
nrow(S5_vs_S1[S5_vs_S1$sig=="significant",]) #only 2 significant genes 
sig_genes_table <- S5_vs_S1[S5_vs_S1$sig=="significant",]
View(sig_genes_table)
S5_vs_S1$gene <- rownames(S5_vs_S1)
colnames(S5_vs_S1)
S5_vs_S1 <- S5_vs_S1[,c("log2FoldChange","pvalue"       ,  "padj","mlogp","sig"      ,      "gene")]
#make master table
master = merge(S5_vs_S1, em_out, by.x=6, by.y=1)
head(master)
write.table(em_out, file="result/master.csv",row.names=FALSE, sep="\t", quote = FALSE)


#visulize #####
#extract top50 genes
top_50 <- master <- master %>%
  dplyr::arrange(padj) %>% dplyr::slice_head(n=50)
em_symbol <- top_50 [,c("gene",sample_info$sample)]

#make hm matrix
#scale 
expr_mat <- as.matrix(em_symbol[, -1])
expr_mat <- t(scale(t(expr_mat)))
rownames(expr_mat) <- em_symbol$gene

#make annotation col 
caf_subtype <- ifelse(
  grepl("CAF-S1", colnames(expr_mat)),
  "S1",
  "S5"
)

annotation_col <- data.frame(
  CAF_Subtype = caf_subtype
)

rownames(annotation_col) <- colnames(expr_mat)


p<-pheatmap::pheatmap(
  expr_mat,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  annotation_col = annotation_col,
  show_colnames = TRUE,
  show_rownames = TRUE
)
png("result/Top50Hm.png", res = 300,width = 10*300, height = 8*300 )
p
dev.off()
