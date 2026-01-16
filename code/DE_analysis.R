library(DESeq2)
library(pheatmap)
library(clusterProfiler)
library(ggplot2)

setwd("D:/Github/CAF_subtypes")



#load data
sc_data<- readRDS("data/integrated_cafs_seurat.rds")
ncol(sc_data) #110 cells

#Pseudo bulk 
cts <- AggregateExpression(sc_data, 
                           group.by = c("subtype", "patient_id"),
                           assays = 'SCT',
                           return.seurat = FALSE)

cts <- cts$SCT

# convert to data.frame
cts <- as.data.frame(cts)
cts = subset(cts,apply(cts, 1, mean) >= 1) #filter low counts 


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
S5_vs_S1 <- S5_vs_S1[,c("log2FoldChange","pvalue"    ,"stat"   ,  "padj","mlogp","sig"      ,      "gene")]
#make master table
master = merge(S5_vs_S1, em_out, by.x=7, by.y=1)
head(master)
write.table(master, file="result/master.csv",row.names=FALSE, sep="\t", quote = FALSE)
write.table(S5_vs_S1, file="result/de.csv",row.names=FALSE, sep="\t", quote = FALSE)


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


#-------------------------------
#Pathway analysis #####

#Load hallmark genes
hallmark <- msigdbr::msigdbr(db_species = "HS",species = "human", collection = "H")
hallmark <- hallmark[,c("gs_name", "gene_symbol")]


#Rank 
x <- S5_vs_S1[!is.na(S5_vs_S1$log2FoldChange) & !is.na(S5_vs_S1$pvalue), , drop = FALSE] #select only non empty rows
x$score <- S5_vs_S1$log2FoldChange
ranks <- x$score
names(ranks) <- x$gene
ranks <- sort(ranks, decreasing = TRUE)
#Run gsea
set.seed(1)
gsea_result <-clusterProfiler::GSEA(
  geneList = ranks,
  TERM2GENE = hallmark,
  pvalueCutoff = 0.25,
  pAdjustMethod = "BH",
  verbose = FALSE
)

#save results 
gsea_result<- gsea_result@result
top_gsea <- gsea_result[gsea_result$p.adjust<0.05,]

write.table(gsea_result,"result/gsea_in_S5.csv",sep="\t", quote = FALSE)
write.table(top_gsea,"result/sig_gsea_in_S5.csv",sep="\t", quote = FALSE)


#Visualize 
#load function####
give_me_pretty_top_gsea <- function(top_gsea, analysis=".",context="."){
  top_gsea$`-logFDR`= -log10(top_gsea$p.adjust)
  p=ggplot(top_gsea, aes(
    x = reorder(ID, NES),
    y = NES,
    fill = NES,
    alpha=`-logFDR`,
    #size = `-logFDR`
  )) +
    geom_point(shape = 21, colour = "black",size=10) +
    scale_fill_gradient2(
      low = "steelblue", mid = "white", high = "red",
      midpoint = 0,
      name = "Enrichment"
    ) + scale_alpha(range = c(0.3, 1),
                    name = "LogFDR") +
    
    
    
    labs(
      x = "Pathway",
      y = "NES",
      title = paste0(analysis,"\n",context)
    ) +
    #my_theme +
    coord_flip()
  return(p)}

p <- give_me_pretty_top_gsea (top_gsea, analysis = "Enriched Pathways in CAF-S5 Vs. CAF-S1", context = "in Late NSCLC")

png("result/topGSEA.png", res=300, width = 300*6, height = 300*6)
p
dev.off()
