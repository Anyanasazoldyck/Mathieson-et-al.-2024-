library(Seurat)
library(stringr)
library(ggplot2)
library(patchwork)
library(sctransform)
library( glmGamPoi )
#set work dir
setwd("C:/Users/testi/Desktop/CAF_subtypes")



sc_list<- readRDS("data/late_nsclc_raw_list.rds")
#merge 
sc_data <- merge(x=sc_list[[1]],y=sc_list[-1])
#sc_data<-SetIdent(sc_data, value = sc_data$patient_id)

#QC####

#Simple  QC eda


sc_data[["percent.mt"]] <- PercentageFeatureSet(sc_data, pattern = "^MT-")
p= VlnPlot(sc_data, features = c("percent.mt"),ncol = 3, pt.size = 0, combine = F)
p2= VlnPlot(sc_data, features = c("nCount_RNA"),ncol = 3, pt.size = 0, combine = F)
p3= VlnPlot(sc_data, features = c("nFeature_RNA"),ncol = 3, pt.size = 0, combine = F)

png("result/Prcent_mt.png", res=300, height = 6*300, width = 10*300)
p 
dev.off()
png("result/nCount_RNA.png", res=300, height = 6*300, width = 12*300)
p2
dev.off()
png("result/nFeature_RNA.png", res=300, height = 6*300, width = 12*300)
p3
dev.off()




plot1 <- FeatureScatter(sc_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
png("result/PreQC_scatter.png", res=300, height = 6*300, width = 15*300)
plot1 + plot2
dev.off()

#filter
#using the threshold advised by Wu et al paper 
sc_data = subset(sc_data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25 & nCount_RNA <25000 )





#normalize using SCT transform ####
sc_data <- SCTransform(sc_data, vars.to.regress = "percent.mt")
#set assay
DefaultAssay(sc_data)<- "SCT"


# Integrate ####
# Run PCA pre integration
sc_data = RunPCA(object = sc_data, assay ="SCT", reduction.name="pca_unintegrated")
ElbowPlot(sc_data, reduction = "pca_unintegrated", ndims = 50)

# decide dims to use, set as a parameter.
dims_to_use = 20

#integration ####
sc_data = harmony::RunHarmony(sc_data,
                     assay.use = "SCT",
                     group.by.vars = c("patient_id"),
                     reduction = "pca_unintegrated", 
                     reduction.save = "pca_integrated")

# re-cluster after integration#####
sc_data = FindNeighbors(sc_data, reduction = "pca_integrated", dims = 1:dims_to_use)
sc_data = FindClusters(sc_data, reduction = "pca_integrated", resolution = 0.1) #As recommended by Wu et al
sc_data = RunUMAP(sc_data, dims = 1:dims_to_use, reduction = "pca_integrated", reduction.name ="umap_integrated")

p<- DimPlot(sc_data, reduction = "umap_integrated", pt.size = 0.5)

png("result/sc_data_umap.png", res = 300, width = 300*6, height = 300*6)
p
dev.off()
#save check point####
#saveRDS(sc_data,"data/integrated_sc_data.rds")

###################CAF IDENTIFICATION #################################
#######################################################################
####  CAF identification and subsetting ####

caf_markers <- c("ITGB1", "PDGFRB", "PDPN", "FAP")   # CD29, PDGFRβ, PDPN, FAP
exclude_marker <- "S100A4"                           # FSP1 
sma_marker <- "ACTA2"                                # αSMA
general_fibroblast_markers <- c("COL1A1", "COL1A2", "DCN")


#Dimplot to identify the Fibroblast cluster 
p <- FeaturePlot(sc_data, features = general_fibroblast_markers, order= T, reduction = "umap_integrated")
p2 <- FeaturePlot(sc_data, features = caf_markers, order= T, reduction = "umap_integrated")
png("result/sc_data_fibroblast.png", res = 300, width = 300*6, height = 300*6)
p
dev.off()
png("result/sc_data_cafs.png", res = 300, width = 300*6, height = 300*6)
p2
dev.off()

#examine the expression of CAF_s1_s2 markerks 

summary (FetchData(sc_data, vars = exclude_marker))
#S100A4      
#Min.   :0.0000  
#1st Qu.:0.0000  
#Median :0.0000  
#Mean   :0.5539  
#3rd Qu.:1.0986  
#Max.   :4.6540  

                               
#ADDmodule score score to explore where CAFs are located

sc_data <- AddModuleScore(sc_data, features = caf_markers, ctrl = 100,
                            name = c("prognostic_cafs"))
sc_data$fibroblast_classification <- ifelse(sc_data$prognostic_cafs1 > 0, "prognostic_cafs","others")
table(sc_data$fibroblast_classification)

p2=DimPlot(sc_data, group.by = "fibroblast_classification",pt.size = 0.5)
png("result/CAFs_ADDModuleScore.png", res = 300, width = 300*6, height = 300*6)
p2
dev.off()

#Subset CAFs ####
fibroblasts_subset <- subset(sc_data, idents = 2)
fibroblasts_subset<- SCTransform(fibroblasts_subset)
# remove FSP1 as neither of these subtype express it ####
summary (FetchData(fibroblasts_subset, vars = exclude_marker))
#S100A4     
#Min.   :0.000  
#1st Qu.:0.000  
#Median :1.099  
#Mean   :1.086  
#3rd Qu.:1.792  
#Max.   :4.554  

fsb1_vlnplot <- VlnPlot(fibroblasts_subset, features = "S100A4")



FSP1_exp <- FetchData(fibroblasts_subset, vars = exclude_marker)
keep_cells <- rownames(FSP1_exp)[FSP1_exp[,1] < 0.5] #ignorable expression 

#how many cells 
cat("fibroblast pre-filtering:", ncol(fibroblasts_subset), "\n")
cat("cells expressing fbs1:", length(keep_cells), "\n")
cat("Percentage:", round(length(keep_cells)/ncol(fibroblasts_subset)*100, 2), "%\n")

#filter cells with < 0.01 expression of FSP1
fibroblasts_subset <- subset(fibroblasts_subset, cells = keep_cells)


#examine prognostic CAF markers
summary (FetchData(fibroblasts_subset, vars = caf_markers))
#ITGB1           PDGFRB           PDPN             FAP        
#Min.   :0.000   Min.   :0.000   Min.   :0.0000   Min.   :0.0000  
#1st Qu.:1.099   1st Qu.:0.000   1st Qu.:0.0000   1st Qu.:0.0000  
#Median :1.609   Median :0.000   Median :0.0000   Median :0.0000  
#Mean   :1.517   Mean   :0.200   Mean   :0.1294   Mean   :0.2252  
#3rd Qu.:2.079   3rd Qu.:0.000   3rd Qu.:0.0000   3rd Qu.:0.0000  
#Max.   :3.689   Max.   :2.197   Max.   :2.3979   Max.   :2.4849  
 

caf_markers_vlnplot <- VlnPlot(fibroblasts_subset, features = caf_markers)
png("result/caf_markers_vlnplot_cluster3.png", res = 300, width = 300*6, height = 300*6)
caf_markers_vlnplot
dev.off()




caf_cells <- WhichCells(
  fibroblasts_subset,
  expression = 
    PDPN > 0.0 &      
    FAP > 0.0    )

cat("Total cells before CAF filtering:", ncol(fibroblasts_subset), "\n")
cat("CAF cells identified:", length(caf_cells), "\n") #110 cell
cat("Percentage of CAFs:", round(length(caf_cells)/ncol(fibroblasts_subset)*100, 2), "%\n") #Percentage of CAFs: 5.68 %

# Subset to CAF cells only
sc_data_cafs <- subset(fibroblasts_subset, cells = caf_cells)

# Save the full dataset before subsetting
saveRDS(sc_data_cafs, "data/sc_data_cafs.rds")

# Continue analysis with CAF subset only
sc_data <- sc_data_cafs
cat("Proceeding with", ncol(sc_data), "CAF cells\n")
#Proceeding with 110 CAF cells


#Calssification and analysis ############################
acta2_vlnplot <- VlnPlot(sc_data, features = "ACTA2")
expr_acta2 <- FetchData(sc_data, vars = "ACTA2")
summary(expr_acta2)
#ACTA2       
#Min.   :0.0000  
#1st Qu.:0.0000  
#Median :0.0000  
#Mean   :0.7537  
#3rd Qu.:1.3144  
#Max.   :3.8067  

#anything in the first and second quantile is probably CAF_s5

Cutoff <-0
sc_data$subtype <- ifelse(
  expr_acta2[,1] > cutoff,
  "CAF_S1",
  "CAF_S5"
)

table(sc_data$subtype)
#CAF_S1 CAF_S5 
#53     57 


# Visualize αSMA expression
p_acta2 <- FeaturePlot(sc_data, features = "ACTA2", order = T, reduction = "umap_integrated", split.by = "subtype")
png("result/acta2_by_subtype_umap.png", res=300, height = 6*300, width = 12*300)
p_acta2
dev.off()

saveRDS(sc_data, "data/integrated_cafs_seurat.rds")


