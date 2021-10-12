## ##CITE-Seq 2 Script####
####Setup####
#Load required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
library(RColorBrewer)
library(writexl)
library(ggridges)
library(clustree)
library(scRepertoire)
library(future)
library(glmGamPoi)

#Alter working capacity
plan()

plan("multiprocess", workers = 4)
options(future.globals.maxSize= 2097152000) # 2Gb

#Initial setup of colour palettes
col = colorRampPalette(brewer.pal(12, 'Set3'))(20)
colbig = colorRampPalette(brewer.pal(12, 'Set3'))(50)

####Load the 10X Cell Ranger output####
#Read the 10x Cell Ranger Output
a_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/GE/A_WT_GE/outs/filtered_feature_bc_matrix")
b_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/GE/B_WT_GE/outs/filtered_feature_bc_matrix")
c_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/GE/C_BCL6_GE/outs/filtered_feature_bc_matrix")
d_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/GE/D_BCL6_GE/outs/filtered_feature_bc_matrix")
f_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/GE/F_E1020K_GE/outs/filtered_feature_bc_matrix")
g_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/GE/G_E1020K_BCL6_GE/outs/filtered_feature_bc_matrix")
h_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/GE/H_E1020K_BCL6_GE/outs/filtered_feature_bc_matrix")


#Add the sample to the cell names, consistent with antibody data below
colnames(a_ge.data)=gsub("-1","_a",colnames(a_ge.data))
colnames(b_ge.data)=gsub("-1","_b",colnames(b_ge.data))
colnames(c_ge.data)=gsub("-1","_c",colnames(c_ge.data))
colnames(d_ge.data)=gsub("-1","_d",colnames(d_ge.data))
colnames(f_ge.data)=gsub("-1","_f",colnames(f_ge.data))
colnames(g_ge.data)=gsub("-1","_g",colnames(g_ge.data))
colnames(h_ge.data)=gsub("-1","_h",colnames(h_ge.data))


#Uppercase the gene names for easier matching later
rownames(a_ge.data)=toupper(rownames(a_ge.data))
rownames(b_ge.data)=toupper(rownames(b_ge.data))
rownames(c_ge.data)=toupper(rownames(c_ge.data))
rownames(d_ge.data)=toupper(rownames(d_ge.data))
rownames(f_ge.data)=toupper(rownames(f_ge.data))
rownames(g_ge.data)=toupper(rownames(g_ge.data))
rownames(h_ge.data)=toupper(rownames(h_ge.data))


head(a_ge.data)
####Load 10X Antibody data####
#Read the 10x Antibody output
a_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/A_WT/umi_count",gene.column=1)
b_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/B_WT/umi_count",gene.column=1)
c_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/C_BCL6/umi_count",gene.column=1)
d_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/D_BCL6/umi_count",gene.column=1)
f_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/F_E1020K_BCL6/umi_count",gene.column=1)
g_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/G_E1020K_BCL6/umi_count",gene.column=1)
h_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/H_E1020K_BCL6/umi_count",gene.column=1)


#Tidy up the rownames from the data
rownames(a_ab.data)=gsub("-[^-]+$","",rownames(a_ab.data),perl=TRUE)
rownames(b_ab.data)=gsub("-[^-]+$","",rownames(b_ab.data),perl=TRUE)
rownames(c_ab.data)=gsub("-[^-]+$","",rownames(c_ab.data),perl=TRUE)
rownames(d_ab.data)=gsub("-[^-]+$","",rownames(d_ab.data),perl=TRUE)
rownames(f_ab.data)=gsub("-[^-]+$","",rownames(f_ab.data),perl=TRUE)
rownames(g_ab.data)=gsub("-[^-]+$","",rownames(g_ab.data),perl=TRUE)
rownames(h_ab.data)=gsub("-[^-]+$","",rownames(h_ab.data),perl=TRUE)

#Add the Sample to the cell names in each sample
colnames(a_ab.data)=paste(colnames(a_ab.data),"_a",sep="")
colnames(b_ab.data)=paste(colnames(b_ab.data),"_b",sep="")
colnames(c_ab.data)=paste(colnames(c_ab.data),"_c",sep="")
colnames(d_ab.data)=paste(colnames(d_ab.data),"_d",sep="")
colnames(f_ab.data)=paste(colnames(f_ab.data),"_f",sep="")
colnames(g_ab.data)=paste(colnames(g_ab.data),"_g",sep="")
colnames(h_ab.data)=paste(colnames(h_ab.data),"_h",sep="")

head(a_ab.data)
####Combine 10X Cell Ranger and Antibody Data into a Suerat Object####
m <- Matrix(nrow = nrow(a_ab.data), ncol = ncol(a_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(a_ab.data)
colnames(m)=colnames(a_ge.data)
common=intersect(colnames(a_ge.data),colnames(a_ab.data))
m[,common]=a_ab.data[,common]
a = CreateSeuratObject(counts = a_ge.data,project="a", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
a[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(b_ab.data), ncol = ncol(b_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(b_ab.data)
colnames(m)=colnames(b_ge.data)
common=intersect(colnames(b_ge.data),colnames(b_ab.data))
m[,common]=b_ab.data[,common]
b = CreateSeuratObject(counts = b_ge.data,project="b", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
b[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(c_ab.data), ncol = ncol(c_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(c_ab.data)
colnames(m)=colnames(c_ge.data)
common=intersect(colnames(c_ge.data),colnames(c_ab.data))
m[,common]=c_ab.data[,common]
c = CreateSeuratObject(counts = c_ge.data,project="c", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
c[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(d_ab.data), ncol = ncol(d_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(d_ab.data)
colnames(m)=colnames(d_ge.data)
common=intersect(colnames(d_ge.data),colnames(d_ab.data))
m[,common]=d_ab.data[,common]
d = CreateSeuratObject(counts = d_ge.data,project="d", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
d[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(f_ab.data), ncol = ncol(f_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(f_ab.data)
colnames(m)=colnames(f_ge.data)
common=intersect(colnames(f_ge.data),colnames(f_ab.data))
m[,common]=f_ab.data[,common]
f = CreateSeuratObject(counts = f_ge.data,project="f", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
f[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(g_ab.data), ncol = ncol(g_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(g_ab.data)
colnames(m)=colnames(g_ge.data)
common=intersect(colnames(g_ge.data),colnames(g_ab.data))
m[,common]=g_ab.data[,common]
g = CreateSeuratObject(counts = g_ge.data,project="g", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
g[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(h_ab.data), ncol = ncol(h_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(h_ab.data)
colnames(m)=colnames(h_ge.data)
common=intersect(colnames(h_ge.data),colnames(h_ab.data))
m[,common]=h_ab.data[,common]
h = CreateSeuratObject(counts = h_ge.data,project="h", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
h[["ADT"]] <- adt_assay

head(a[[]])
#####Process samples as one####
experiments=c(a,b,c,d,f,g,h)
experiment_names=c("a","b","c","d","f","g","h")

experiment<-merge(x= a, y=c(b,c,d,f,g,h))

experiment
str(experiment)
head(experiment[[]])
####Quality control, filtering, normalisation and scaling####
#Mitochondrial QC metrics
experiment[["percent.mt"]] <- PercentageFeatureSet(experiment, pattern = "^MT-")

#Remove where nCount_ADT = 0
DefaultAssay(experiment) <- "ADT"
experiment <- subset(experiment, nCount_ADT >0 )
DefaultAssay(experiment) <- "RNA"

#Visualize QC metrics as violin plot
RNA_QC <- VlnPlot(experiment, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ADT_QC <- VlnPlot(experiment, features = c("nFeature_ADT", "nCount_ADT", "percent.mt"))

#FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 = FeatureScatter(experiment, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()  +
  ylab("% of mitochondrial genes") +
  xlab("UMI counts") + 
  geom_hline(yintercept = 5) 

plot2 = FeatureScatter(experiment, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() +
  ylab("Number of genes") +
  xlab("UMI counts") + 
  geom_hline(yintercept = 200) 

plot1 + plot2

#Generally aim to filter out unique feature counts over 2,500 and less than 200; and percent.mt over 5%
filter_seurat = function(seurat_object){
  
  message("Performing filter by number of genes and mitochondrial percentage.")
  seurat_object = subset(seurat_object, subset = nFeature_RNA > 200  & percent.mt < 5 & nFeature_RNA < 2500)
  message("Now the object has ", dim(seurat_object)[1], " genes and ", dim(seurat_object)[2], " cells.")
  
  
  return(seurat_object)
}

experiment <-  filter_seurat(experiment)

#Normalise dataset - default
experiment <- NormalizeData(experiment, normalization.method = "LogNormalize", scale.factor = 10000)

#Normalise dataset - SCTransform
experiment = SCTransform(experiment, verbose = TRUE)

#Find variable features
experiment <- FindVariableFeatures(experiment, selection.method = "vst")
?FindVariableFeatures
top20 = head(VariableFeatures(experiment), 20)

plot3 = VariableFeaturePlot(experiment)
plot4 = LabelPoints(plot = plot3, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)
plot4

#Scale data
experiment <- ScaleData(experiment)

####Dimensionality reduction -PCA####
#Perform linear dimensional reduction (PCA)
experiment <- RunPCA(experiment, verbose = FALSE, features = VariableFeatures(object = experiment))

##Visualise PCA results
print(experiment[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(experiment, dims = 1:4, reduction = "pca", nfeatures = 15)

plot5 <- DimPlot(experiment, reduction = "pca", dims = c(1,2))
plot6 <- DimPlot(experiment, reduction = "pca", dims = c(1,3))
plot5 + plot6

DimHeatmap(experiment, dims = 1:6, cells = 500, balanced = TRUE)

####Determine dimensionality of the dataset - How many principal components should be included to capture the majority of variance?####
pca_variance <- experiment@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #23

####Determine number of clusters####
##luster Tree Analysis
clustree(experiment, prefix = "RNA_snn_res.") +
  theme(legend.position="bottom")

#Cluster the cells
experiment <- FindNeighbors(experiment, dims = 1:30)
experiment <- FindClusters(experiment, resolution = 1.0, verbose = FALSE)
experiment <- RunUMAP(experiment, dims = 1:30)
DimPlot(experiment, label = TRUE, cols=colbig) +  ggtitle("RNA Clustering")

####Scale antibody data####
DefaultAssay(experiment) <- "ADT"
VariableFeatures(experiment) <- rownames(experiment[["ADT"]])
experiment <- NormalizeData(experiment, normalization.method = "CLR", margin = 2)
experiment <- ScaleData(experiment)
experiment <- RunPCA(experiment,reduction.name = 'apca')

#Visualise antobody PCA
print(experiment[["apca"]], dims = 1:10, nfeatures = 5)
Plot_13 <- VizDimLoadings(experiment, dims = 1:4, reduction = "apca", nfeatures = 15)
Plot_13

####Combine into wnn plot####
experiment <- FindMultiModalNeighbors(
  experiment, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)
#Determing resolution (cluster number) for wnn plot?
clustree(experiment, prefix = "wsnn_res.") +
  theme(legend.position="bottom")


#UMPA plots for RNA, ADT and WNN
experiment <- RunUMAP(experiment, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                      reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
experiment<- RunUMAP(experiment, reduction = 'apca', dims = 1:18, assay = 'ADT', 
                     reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
experiment <- RunUMAP(experiment, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
experiment <- FindClusters(experiment, graph.name = "wsnn", algorithm = 3, resolution = 1.0, verbose = TRUE)


DefaultAssay(experiment) <- "RNA"

p1=DimPlot(experiment, label = TRUE,cols=colbig,reduction = "rna.umap", label.size = 2.5) + NoLegend()
p2=DimPlot(experiment, label = TRUE,cols=colbig,reduction = "adt.umap", label.size = 2.5) + NoLegend()
p3=DimPlot(experiment, label = TRUE,cols=colbig, reduction = "wnn.umap", label.size = 2.5) + NoLegend()


p1|p2|p3

p1
p2
p3

###Umap-wnn by mouse
plot_mouse <- DimPlot(experiment, label = TRUE,reduction = "wnn.umap", label.size = 2.5, group.by = "orig.ident")
plot_mouse

###Match the RNA Names to the Antibodies, this should be checked
list1=c(rownames(c1_ab.data))
list2=c("PTPRC","FAS","CD19","IGHM","CR2","FCER2A","CD93","CD83","CD86","IGHD","CD8A","SELL","CD44","CD4","CXCR5","PDCD1","IL2RA","CD274","PDCD1LG2","CTLA4","CD80","CD40","CD69","ICOS","CD38","TNFRSF18")

####Analysis of clusters####
##Different plotting options
DefaultAssay(experiment) <- "RNA"
FeaturePlot(experiment, features = c("CD19", "CD4", "CD8A", "PRDM1", "PPBP", "NKG7", "CST3", "FOXP3", "B220"), reduction = "wnn.umap")

RidgePlot(experiment, features = c("CD19", "CYP11A1"), ncol = 2)

FeaturePlot(experiment, features = "CD4", reduction = "wnn.umap")

##Finding all the markers
experiment.markers <- FindAllMarkers(experiment, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
experiment.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(experiment, features = top10$gene) + NoLegend()

##DE genes of individual clusters
Cluster_2 <- FindMarkers(experiment, ident.1 = 2, assay = "RNA")
Cluster_2_adt <- FindMarkers(experiment, ident.1 = 2, assay = "ADT")

Cluster_12 <- FindMarkers(experiment, ident.1 = 12, assay = "RNA")
Cluster_12_adt <- FindMarkers(experiment, ident.1 = 12, assay = "ADT")

##Subsetting unknown cluster
Unknown_cells <- subset(experiment, idents = c(2, 12, 27, 33))
Unknown_cells <- FindClusters(Unknown_cells, resolution = 0.8, verbose = FALSE, graph.name = "wsnn")
Unknown_cells <- RunUMAP(Unknown_cells, dims = 1:30, reduction.name = "unknown.umap")
DimPlot(Unknown_cells, label = TRUE, cols=colbig, reduction = "unknown.umap", label.size = 2.5) + NoLegend()
FeaturePlot(Unknown_cells, "SOX4", reduction = "unknown.umap")

###Clonotype analysis
#Data preparation
contig_list <- list(c1_cl.data, c2_cl.data, d1_cl.data, d2_cl.data)
head(contig_list[[1]])

combined <- combineBCR(contig_list, samples = c("MouseC1", "MouseC2", "MouseD1", "MouseD2"), ID = c("c1", "c2", "d1", "d2"))
str(combined)
head(combined[[1]])

##Data visualisation
#Percent/total number of unique clonotypes 
quantContig(combined, cloneCall = "gene+nt", scale = T) #percent of unique clonotypes of total size of the size of clonotyeps
quantContig(combined, cloneCall = "gene+nt", scale = F) #number of uniqe clonotypes

quantContig(combined, cloneCall = "gene+nt", scale = T, chain = "IGL")#chain argument not working

#Abundance of clonotypes
Abundance_clonotypes <- abundanceContig(combined, cloneCall = "gene", scale = F, exportTable = T)
Abundance_clonotypes <- Abundance_clonotypes %>%
  arrange(desc(Abundance))
Abundance_clonotypes

abundanceContig(combined, cloneCall = "gene", scale = T)

#Length of clonotypes
lengthContig(combined, cloneCall = "aa")
lengthContig(combined, cloneCall = "nt")

#Compare clonotypes
compareClonotypes(combined, samples = c("MouseC1", "MouseD1"), cloneCall = "aa", graph = "alluvial") #Does not work

#Visualise Gene Usage
