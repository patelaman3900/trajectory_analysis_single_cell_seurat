#---------------------------------------------sc Trajectory Analysis-------------------------------------------------------------------

#WHAT is trajectory analysis?
#To check the lineage of cell differentitation by using the expression of gene and its changes
#To place the cell in the given pathway of development

#What is pseudotime?
#Cell at earlier stage (smaller pseudotime);Cell at later stage (Higher pseudotime)
#Abstract unit of progress
#Measure of progress an individual cell has made through a process such as cell differentiation

#WHEN to perform trajectory analysis?
#General assumptions:
#Biological process of interest is dynamic and its mechanism
#Prior knowledge and/or evidence that a trajectory exists
#Appropriate cells are sampled
#Data sampled to sufficient depth ensuring the presence of continuum of states among cells

#WHICH trajectory inference method to choose?
#Disconnected trajectories?
#Prior knowledge?
#Particular topology expected?
#User Friendly?
#Compatibilities
#Performance
#dynverse()- collection of R packages aimed at supporting the trajectory inference(TI) community on multiple levels

#OVERALL WORKFLOW
#ScRNA-seq dataset
#Pre-process data (Normalize , Remove batch effects)
#Non-linear dimensionality reduction (1-SNE, UMAP)
#Cluster cells
#Compare clusters (Identify top markers, Targeted contrasts)
#Trajectory analysis

#monocle3 WORKFLOW
#preprocess_cds()
#reduce_dimension()
#cluster_cells()
#top_markers()
#learn_graph()

#Seurat+monocle3 WORKFLOW
#Seurat workflow steps + findMarkers()/automatic annotation tools

#NOTE:-
#cell_data_set object structure
#monocle3 requires cell_data_set class
#cell_data_set class is derived from SingleCellExperiment Object
#Takes the count matrix produced by cell ranger pipelines(10X genomics)
#cell_data_set require 3 input files
#1. Expression Matrix(Rows are genes,columns are cells)       
#2. Cell metadata dataframe (Rows are cells, colums are attributes(cell types,culture conditions,day capture))  
#3. Gene metadat dataframe (Rows are genes, columns are attribute of the gene(GC content,etc))

#DATA Detail
#7,551 human blood cells were profiled using scRNA-seq (STRT-seq), covering 32 immunophenotypic cell types from 21 healthy donors.
#For the demo today, we will be subsetting B cells and progenitors (1,448 cells).

#Goal of the analysis:
#Construct a trajectory
#Order cells in pseudotime
#Find genes that change expression as cells progress along a trajectory

#HOW to perform this analysis?

#Required R Packages
#monocle3
#Seurat
#SeuratWrappers
#tidyverse

#SCRIPT

# script to perform trajectory analysis
# https://www.nature.com/articles/s41467-019-10291-0
#setpath 
setwd("D:/AMAN_PhD/Script/scRNA_seq/Seurat/scTrajectory_Analysis")

set.seed(1234)

#install packages
install.packages("BiocManager")
BiocManager::install(version = "3.19")

install.packages(c("Seurat", "ggplot2", "tidyverse", "remotes"))
install.packages(c("R.utils", "igraph", "uwot", "FNN", "proxy"))
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/grr/grr_0.9.5.tar.gz",
  repos = NULL,
  type = "source"
)


remotes::install_github(
  "cole-trapnell-lab/monocle3",
  ref = "master",
  dependencies = TRUE
)


install.packages("remotes")
remotes::install_github("satijalab/seurat-wrappers")



#load libraries
library(grr)
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)


# read in data
markers <- read.delim('D:/AMAN_PhD/Script/scRNA_seq/Seurat/scTrajectory_Analysis/ABC_Marker.txt', header = T) # gene metadata
metadata <- read.delim('D:/AMAN_PhD/Script/scRNA_seq/Seurat/scTrajectory_Analysis/ABC_Meta.txt', header = T) # cell metadata
expr <- read.delim('D:/AMAN_PhD/Script/scRNA_seq/Seurat/scTrajectory_Analysis/ABC_umi_matrix_7551_cells.csv', header = T, sep = ',') # expression matrix



# create seurat object ---------------
expr.t <- t(expr)    #transpose the expression matrix
seu.obj <- CreateSeuratObject(counts = expr.t)  
View(seu.obj@meta.data)
seu.obj@meta.data <- merge(seu.obj@meta.data, metadata, by.x = 'row.names', by.y = 'cell_id')
View(seu.obj@meta.data)
seu.obj@meta.data <- seu.obj@meta.data %>% 
  column_to_rownames(var = 'Row.names')

Idents(seu.obj) <- factor(Idents(seu.obj), levels = levels(Idents(seu.obj)))
names(Idents(seu.obj)) <- colnames(seu.obj)


cells <- colnames(seu.obj)

counts <- GetAssayData(
  seu.obj,
  assay = "RNA",
  layer = "counts"
)

# reorder explicitly
counts <- counts[, cells]

seu.obj <- SetAssayData(
  seu.obj,
  assay = "RNA",
  layer = "counts",
  new.data = counts
)

validObject(seu.obj)





seu.obj$mitopercent <- PercentageFeatureSet(seu.obj, pattern = '^MT-')
view(seu.obj$mitopercent)
seu.obj.filtered <- subset(seu.obj, subset = nCount_RNA > 800 &
                             nFeature_RNA > 500 &
                             mitopercent < 10)


# subset my seurat object - B cells

unique(seu.obj.filtered@meta.data$population)

Idents(seu.obj.filtered) <- seu.obj.filtered$population
b.seu <- subset(seu.obj.filtered, idents = "b")
b.seu
unique(b.seu@meta.data$redefined_cluster)

# pre-processing using seurat
b.seu <- NormalizeData(b.seu)
b.seu <- FindVariableFeatures(b.seu)
b.seu <- ScaleData(b.seu)
b.seu <- RunPCA(b.seu)
b.seu <- FindNeighbors(b.seu, dims = 1:30)
b.seu <- FindClusters(b.seu, resolution = 0.9)
b.seu <- RunUMAP(b.seu, dims = 1:30, n.neighbors = 50)

a1 <- DimPlot(b.seu, reduction = 'umap', group.by = 'redefined_cluster', label = T)
a2 <- DimPlot(b.seu, reduction = 'umap', group.by = 'seurat_clusters', label = T)
a1|a2


# MONOCLE3 WORKFLOW ---------------------
# monocle3 requires cell_data_set object
# convert seurat object to cell_data_set object for monocle3





# ...1 Convert to cell_data_set object ------------------------

cds <- as.cell_data_set(b.seu)
cds

cds <- cluster_cells(cds)

# to get cell metadata
colData(cds)
# to gene metdata
fData(cds)
rownames(fData(cds))[1:10]

# since it misses the gene_short_name column, let's add it
fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)



# ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------
# let's use the clustering information have

# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)


cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 

list_cluster <- b.seu@active.ident
cds@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <- b.seu@reductions$umap@cell.embeddings



# plot

cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
                            color_cells_by = "redefined_cluster",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan')) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names



# ...3. Learn trajectory graph ------------------------
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = 'redefined_cluster',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)


# ...4. Order the cells in pseudotime -------------------

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 5]))

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)

# cells ordered by monocle3 pseudotime

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(redefined_cluster, monocle3_pseudotime, median), fill = redefined_cluster)) +
  geom_boxplot()




# ...5. Finding genes that change as a function of pseudotime --------------------
deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()

FeaturePlot(b.seu, features = c('E2F2', 'STMN1', 'CD52'))


# visualizing pseudotime in seurat

b.seu$pseudotime <- pseudotime(cds)
Idents(b.seu) <- b.seu$redefined_cluster
FeaturePlot(b.seu, features = "pseudotime", label = T) 



#DATA SOURCE

#Link to code:
#https://github.com/kpatel427/YouTubeTutorials/blob/main/TI_monocle3.R

#Data for analysis:
#http://scrna.sklehabc.com/

#Alternate Data Link:
#https://drive.google.com/file/d/1CJ9VSrUCoqPsUI1jrdm2nrLRawI04xZ1/view

#Publication associated with the data:
#https://doi.org/10.1093/nsr/nwaa180

#Monocle3 tutorial:
#https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/

#R package collection for Trajectory Inference:
#https://dynverse.org/

#Publication comparing various Trajectory Inference methods:
#https://www.biorxiv.org/content/10.1101/276907v1.full.pdf



#---------------------------------------------sc Trajectory Analysis-------------------------------------------------------------------
