################################################################################
#  Myeloid cells in astrocytoma and oligodendroglioma - analysis script        #
################################################################################
#
# OVERVIEW. This script characterises the myeloid compartment (macrophages and
# microglia) of glioma tumours using two single-cell RNA-seq cohorts,
# astrocytoma and oligodendroglioma. Starting from the cells labelled
# "Macrophage" in the global immune annotation of each cohort, it cleans and
# integrates each cohort, merges them, defines four myeloid clusters, and then
# produces the panels of Figure 6 and of its supplementary figures (marker-gene
# heatmaps, cluster proportions, UMAPs, transcription-factor activity and
# Hallmark pathway enrichment).
#
# WORKFLOW. The sections below are meant to be run in this order:
#
# 1. Pre-processing: per-cohort normalisation, integration across samples and
#    removal of contaminating (neuronal / T-cell) clusters.
# 2. Merged dataset and metadata
# 3. Program scoring: UCell scores of the myeloid programs described by Miller
#    et al.; each cell is labelled with its top-scoring program.
# 4. Cluster definition: removal of one outlier sample, re-clustering and
#    annotation of the four myeloid clusters (Microglia, Inflammatory,
#    Macrophage, SPP1+ Macrophage).
# 5. Heatmaps of cluster marker genes (supplementary, then main).
# 6. Proportion plots of the clusters (Fig. 6C, 6E, 6F).
# 7. UMAP panels (Fig. 6A main, 6B supp).
# 8. Transcription-factor analysis from SCENIC+ (Fig. 6D, 6E supp).
# 9. Hallmark over-representation analysis (Fig. 6A supp).
#
# INPUTS. astrocytoma_all_immune.rds and oligodendroglioma_all_immune.rds
# (Seurat objects of all immune cells; the column cell.type.sr.cl.imm.postmerge
# is used to keep the macrophages); Myeloid_NMF_Average_Gene_Spectra.txt (gene
# weights of the myeloid NMF programs of Miller et al.);
# direct_gene_based_AUC_astro.tsv.gz and direct_gene_based_AUC_oligo.tsv.gz
# (per-cell SCENIC+ regulon AUC scores, used in the TF section).
#
# OUTPUTS.
# all.myeloid.filt.rds (final object with the cluster annotation);
# cell_annotations_minus_S0009.tsv (cluster label and UMAP coordinates of every
# cell, re-used by the proportion and UMAP sections); heatmap_*_to_plot.tsv
# (cells and genes shown in the heatmaps); the figures, written in the working
# directory, in output/, TF_heatmaps/ and results_myeloid_ora/.
#
# PRACTICAL NOTES. Sample names follow the pattern
# <disease>_<patient>_<status>_<sample>, e.g. oligo_P0004_ini_S0009 (status: ini
# = initial, rec = recurrence); this pattern is parsed in several places. The
# folder output/ must exist before running. rm(list = ls()) is called twice
# to free memory, so the following sections reload what they need from disk.
#
################################################################################

# Packages
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(UCell)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(patchwork)
library(msigdbr)
library(viridis)
library(lme4)
library(lmerTest)
library(harmony)      
library(batchelor)    
library(glmGamPoi)    
library(presto)       
pal = viridis(n=50)



# Dataset
astro.ss= readRDS('astrocytoma_all_immune.rds')   # All myeloid cells as define as Macrophage in the other part of the study in astrocytoma
oligo.ss= readRDS('oligodendroglioma_all_immune.rds')   # All myeloid cells as define as Macrophage in the other part of the study in oligodendroglioma

astro.ss <- subset(astro.ss, cells = colnames(astro.ss)[astro.ss$cell.type.sr.cl.imm.postmerge == "Macrophage"])
oligo.ss <- subset(oligo.ss, cells = colnames(oligo.ss)[oligo.ss$cell.type.sr.cl.imm.postmerge == "Macrophage"])

program_matrix = read.table('Myeloid_NMF_Average_Gene_Spectra.txt') # Program matrice from Miller et al.


###########################
#Pre-processing functions #
###########################

process_myeloid_integration <- function(obj, sample_col = "orig.ident",
                                        exclude_samples = NULL,
                                        rpca_k_weight,
                                        res_low = 0.3, res_high = 2.0) {

  obj[["RNA"]] <- split(obj[["RNA"]], f = obj[[sample_col, drop = TRUE]])
  obj <- SCTransform(obj, assay = "RNA", vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = FALSE)

  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj, dims = 1:50, reduction = "pca", graph.name = "unintegrated_nn")
  obj <- FindClusters(obj, resolution = res_low, cluster.name = paste0("unintegrated_clusters", res_low), graph = "unintegrated_nn")
  obj <- RunUMAP(obj, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")

  obj <- IntegrateLayers(object = obj, method = HarmonyIntegration,
                         orig.reduction = "pca", new.reduction = "integrated.har",
                         verbose = FALSE, assay = "RNA")
  obj <- FindNeighbors(obj, reduction = "integrated.har", dims = 1:50, graph.name = "har_nn")
  obj <- FindClusters(obj, resolution = res_high, cluster.name = paste0("har_clusters", res_high), graph = "har_nn")
  obj <- RunUMAP(obj, reduction = "integrated.har", dims = 1:50, reduction.name = "umap.har")

  obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,
                         new.reduction = "integrated.mnn", verbose = FALSE, batch = obj[[sample_col, drop = TRUE]])
  obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:50, graph.name = "fmnn_nn")
  obj <- FindClusters(obj, resolution = res_low, cluster.name = paste0("fmnn_clusters", res_low), graph = "fmnn_nn")
  obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.fmnn")

  if (!is.null(exclude_samples)) {
    keep_cells <- colnames(obj)[!(obj[[sample_col, drop = TRUE]] %in% exclude_samples)]
    obj <- subset(obj, cells = keep_cells)
  }

  DefaultAssay(obj) <- "RNA"
  obj[["SCT"]] <- NULL
  obj <- SCTransform(obj, assay = "RNA", vars.to.regress = c("percent.mt", "nCount_RNA"))
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- IntegrateLayers(object = obj, method = RPCAIntegration, new.reduction = "integrated.rpca",
                         normalization.method = "SCT", verbose = FALSE, k.weight = rpca_k_weight)
  obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:50, graph.name = "rpca_nn")
  obj <- FindClusters(obj, resolution = res_low, cluster.name = paste0("rpca_clusters", res_low), graph = "rpca_nn")
  obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:50, reduction.name = "umap.rpca")

  obj
}


remove_clusters <- function(obj, cluster_col, clusters_to_remove) {
  Idents(obj) <- cluster_col
  cells_to_keep <- WhichCells(obj, idents = clusters_to_remove, invert = TRUE)
  subset(obj, cells = cells_to_keep)
}

identify_clusters <- function(obj, cluster_col, reduction,
                              marker_genes = c("RBFOX3", "SNAP25", "STMN2", "CD3E", "CD247")) {
  Idents(obj) <- cluster_col
  print(DimPlot(obj, reduction = reduction, label = TRUE))
  print(FeaturePlot(obj, features = marker_genes, reduction = reduction))
  obj <- PrepSCTFindMarkers(obj, assay = "SCT")
  FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25)
}

plot_before_after <- function(obj_before, obj_after, cluster_col, reduction, label) {
  print(DimPlot(obj_before, reduction = reduction, group.by = cluster_col, label = TRUE) +
          ggtitle(paste(label, "- before")))
  print(DimPlot(obj_after, reduction = reduction, group.by = cluster_col, label = TRUE) +
          ggtitle(paste(label, "- after")))
}

#######################
# Astrocytoma dataset #
#######################

astro.ss <- process_myeloid_integration(astro.ss, rpca_k_weight = 79)
astro_markers <- identify_clusters(astro.ss, "har_clusters2", "umap.har")

clusters_check <- c("15","16")
astro_markers %>%
  filter(cluster %in% clusters_check, p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  print(n = Inf)

astro_macro_sub <- remove_clusters(astro.ss, "har_clusters2", c("15","16")) # neuron cluster 
plot_before_after(astro.ss, astro_macro_sub, "har_clusters2", "umap.har", "Astro")

#############################
# Oligodendroglioma dataset #
#############################

oligo.ss <- process_myeloid_integration(
  oligo.ss,
  exclude_samples = c("oligo_P0006_rec1_S0014", "oligo_P0024_rec1_S0047", "oligo_P0001_rec2_S0003"), #low number of cells
  rpca_k_weight = 56
)
oligo_markers <- identify_clusters(oligo.ss, "har_clusters2", "umap.har")

clusters_check <- c("20","16","22")

oligo_markers %>%
  filter(cluster %in% clusters_check, p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  print(n = Inf)

oligo_macro_sub <- remove_clusters(oligo.ss, "har_clusters2", c("16","20")) # neuron and lymphocyte T clusters 
plot_before_after(oligo.ss, oligo_macro_sub, "har_clusters2", "umap.har", "Oligo")

##########################################
# Merged dataset                         #
##########################################

merged.ss <- merge(astro_macro_sub, oligo_macro_sub)
DefaultAssay(merged.ss) <- "SCT"

merged.ss <- NormalizeData(merged.ss)
merged.ss <- FindVariableFeatures(merged.ss)
merged.ss <- ScaleData(merged.ss)
merged.ss <- RunPCA(merged.ss)

merged.ss <- IntegrateLayers(object = merged.ss, method = HarmonyIntegration,
                              orig.reduction = "pca", new.reduction = "integrated.har",
                              verbose = FALSE, assay = "RNA")
merged.ss <- FindNeighbors(merged.ss, reduction = "integrated.har", dims = 1:50, graph.name = "har_nn")
merged.ss <- RunUMAP(merged.ss, reduction = "integrated.har", dims = 1:50, reduction.name = "umap.har")
merged.ss <- FindClusters(merged.ss, graph.name = "har_nn", resolution = 0.2, cluster.name = "clusters_0.2")


########################################
# Add some metadata                    #
########################################
merged.ss$disease <- ifelse(
  grepl("astro", merged.ss$orig.ident, ignore.case = TRUE),
  "astro",
  ifelse(
    grepl("oligo", merged.ss$orig.ident, ignore.case = TRUE),
    "oligo",
    NA
  )
)

parts <- strsplit(as.character(merged.ss$orig.ident), "_")

merged.ss$patient <- vapply(parts, function(x) x[2], character(1))
merged.ss$status <- vapply(parts, function(x) x[3], character(1))
merged.ss$sample <- vapply(parts, function(x) x[4], character(1))

#DimPlot(merged.ss, reduction = "umap.har", group.by = "disease", label = TRUE)

##################################################
#  Scoring with UCell the Miller et al. programs #
##################################################
gene_sets <- list()

# Use the top genes
for (prog in colnames(program_matrix)) {
  
  scores = program_matrix[, prog]
  names(scores) = rownames(program_matrix)
  
  top_genes = names(sort(scores, decreasing = TRUE))[1:100]
  top_genes = top_genes[top_genes %in% rownames(merged.ss)]
  
  if (length(top_genes) >= 10) {
    gene_sets[[prog]] = top_genes
  }
}

#Scoring using UCell
all.myeloid <- AddModuleScore_UCell(
  merged.ss,
  features = gene_sets,
  assay = 'SCT',
  name = "Program"
)

# Plot
FeaturePlot(
  all.myeloid,
  reduction = "umap.har",
  features = grep("Program$", colnames(all.myeloid@meta.data), value = TRUE)
) & scale_colour_gradientn(colours = pal, na.value = "grey90")


# Change the names
cols <- grep("Program$", colnames(all.myeloid@meta.data), value = TRUE)
colnames(all.myeloid@meta.data)[
  match(cols, colnames(all.myeloid@meta.data))
] <- sub("Program$", "", cols)


# Add a label
program_cols <- c("Microglia","Complement_Immunosuppressive","Neutrophil",
                  "Scavenger_Immunosuppressive","G2_M","IFN_Response",
                  "Inflammatory_microglia","Macrophage",
                  "IL1B_CD83_TNF_Inflammatory","cDC","G1_S","Hypoxia",
                  "HS_UPR","Monocyte")
all.myeloid$label_program <- apply(
  all.myeloid@meta.data[, program_cols, drop = FALSE],
  1,
  function(x) names(x)[which.max(x)]
)

#saveRDS(all.myeloid, "all.myeloid.rds")

############################################################
# Data sample removing and attribution of the cluster name #
############################################################
#all.myeloid = readRDS("all.myeloid.rds")

N.per.cat = 20
Ncell = 1000

#remove strange S0009 oligo samples which has high interferon response
samps.to.rem = c( "oligo_P0004_ini_S0009")

cells.to.keep <- colnames(all.myeloid)[
  !all.myeloid$orig.ident %in% samps.to.rem
]

all.myeloid.filt <- subset(all.myeloid, cells = cells.to.keep)

all.myeloid.filt <- IntegrateLayers(
  object = all.myeloid.filt, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.har",
  verbose = FALSE, assay = "RNA"
)
all.myeloid.filt <- FindNeighbors(all.myeloid.filt, reduction = "integrated.har", dims = 1:30, graph.name = "har_nn")
all.myeloid.filt <- RunUMAP(all.myeloid.filt, reduction = "integrated.har", dims = 1:30, reduction.name = "umap.har")
res = 0.3
all.myeloid.filt <- FindClusters(all.myeloid.filt, resolution = res, cluster.name = paste0("har_clusters",res), graph = "har_nn", algorithm = 1)
DimPlot(all.myeloid.filt, group.by=c(paste0("har_clusters",res),"disease","patient","grade","orig.ident"), 
        reduction = "umap.har", shuffle = T)


FeaturePlot(all.myeloid.filt, features = c("Microglia","Complement_Immunosuppressive","Neutrophil","Scavenger_Immunosuppressive","G2_M","IFN_Response",
                                           "Inflammatory_microglia","Macrophage","IL1B_CD83_TNF_Inflammatory","cDC","G1_S","Hypoxia","HS_UPR","Monocyte","SPP1"),
            reduction = "umap.har")


all.myeloid.filt = PrepSCTFindMarkers(all.myeloid.filt)
Idents(all.myeloid.filt) <- paste0("har_clusters",res)
markers = FindAllMarkers(all.myeloid.filt)

genes.to.plot = unlist(lapply(c(0:3), function(x) {
  if(x==3){
    markers[markers$cluster==x & markers$avg_log2FC>0,"gene"][1:(N.per.cat+1)]
  }else{
    markers[markers$cluster==x & markers$avg_log2FC>0,"gene"][1:N.per.cat]
  }
}))
#remove second CD163 since it is a marker of clusters 3 and 4
inds = which(genes.to.plot=="CD163")
genes.to.plot = genes.to.plot[-inds[2]]

#subsample 1000 cells per cluster to plot, plus all cycling cells
cells.to.plot.list = lapply(c(0:3), function(x) {
  all.cells.cl = names(all.myeloid.filt$har_clusters0.3[which(all.myeloid.filt$har_clusters0.3==x)])
  all.cells.cl[sample(length(all.cells.cl),min(Ncell,length(all.cells.cl)))]
})
cells.to.plot = unlist(cells.to.plot.list)

#write lists of cells and genes for figure
write.table(cells.to.plot, "heatmap_main_cells_to_plot.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(genes.to.plot, "heatmap_main_genes_to_plot.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#write full list of cells, umap coordinates, and cluster assignments
all.myeloid.filt$cluster_annotation = all.myeloid.filt$har_clusters0.3
all.myeloid.filt$cluster_annotation = as.numeric(as.character(all.myeloid.filt$cluster_annotation))
all.myeloid.filt$cluster_annotation[all.myeloid.filt$cluster_annotation==0] = "Microglia"
all.myeloid.filt$cluster_annotation[all.myeloid.filt$cluster_annotation==1] = "Inflammatory"
all.myeloid.filt$cluster_annotation[all.myeloid.filt$cluster_annotation==2] = "Macrophage"
all.myeloid.filt$cluster_annotation[all.myeloid.filt$cluster_annotation==3] = "SPP1+ Macrophage"
df.to.write = cbind(data.frame(cell.name = colnames(all.myeloid.filt),
                               cluster = all.myeloid.filt$cluster_annotation),
                    Embeddings(all.myeloid.filt, reduction = "umap.har"))

write.table(df.to.write, "cell_annotations_minus_S0009.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
colnames(all.myeloid@meta.data)[colnames(all.myeloid@meta.data) == "cluster_annotation"] <- "cluster"

saveRDS(all.myeloid.filt, "all.myeloid.filt.rds")
#########################
# Data for heatmap supp #
#########################
cells.main <- read.table('heatmap_main_cells_to_plot.tsv', header = FALSE, sep = "\t")[,1]
genes.main <- read.table('heatmap_main_genes_to_plot.tsv', header = FALSE, sep = "\t")[,1]

set.seed(42)
cells.s0009 <- colnames(all.myeloid)[all.myeloid$orig.ident %in% samps.to.rem]
cells.ifn <- cells.s0009[sample(length(cells.s0009), min(Ncell, length(cells.s0009)))]

cells.secondary <- c(cells.main, cells.ifn)
labels.secondary <- c(rep(c("Microglia", "Inflammatory", "Macrophage", "SPP1+ Macrophage"), each = Ncell),
                      rep("Interferon-response", length(cells.ifn)))

all.myeloid <- PrepSCTFindMarkers(all.myeloid)
Idents(all.myeloid) <- ifelse(all.myeloid$orig.ident %in% samps.to.rem, "S0009", "other")
mk.s0009 <- FindMarkers(all.myeloid, ident.1 = "S0009", ident.2 = "other", only.pos = TRUE)

genes.main <- unique(genes.main)
genes.ifn  <- head(setdiff(rownames(mk.s0009), genes.main), 100 - length(genes.main))
genes.secondary <- c(genes.main, genes.ifn)

stopifnot(length(genes.secondary) == 100, !any(duplicated(genes.secondary)))

write.table(data.frame(cell = cells.secondary, cluster = labels.secondary),
            "heatmap_secondary_cells_to_plot.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(genes.secondary, "heatmap_secondary_genes_to_plot.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

rm(list = ls())
gc()

################################
# Heatmap gene expression Supp #
################################

all.myeloid <- readRDS('all.myeloid.rds')
cells.to.plot <- read.table('heatmap_secondary_cells_to_plot.tsv', header = TRUE, sep = "\t")[,1]
genes.to.plot <- read.table('heatmap_secondary_genes_to_plot.tsv', header = FALSE, sep = "\t")[,1]

all.genes.plot <- c(genes.to.plot)

Ncell = 1000
N.per.cat = 20

mat <- all.myeloid@assays$SCT@data
mat.scaled <- t(scale(t(mat)))

col.annot.df <- data.frame(
  cluster = c(
    rep("Microglia", Ncell),
    rep("Inflammatory", Ncell),
    rep("Macrophage", Ncell),
    rep("SPP1+ Macrophage", Ncell),
    rep("Interferon-response", Ncell)
  )
)
rownames(col.annot.df) <- cells.to.plot

col.annot.df$cluster <- factor(col.annot.df$cluster,
                               levels = c("Microglia", "Inflammatory", "Macrophage","SPP1+ Macrophage","Interferon-response"))

meta.sub <- all.myeloid@meta.data[rownames(col.annot.df), ]
meta.sub$disease <- trimws(as.character(meta.sub$disease))

col.annot.df$grade  <- factor(as.character(meta.sub$grade), levels = c("2", "3", "4"))
col.annot.df$cancer <- factor(meta.sub$disease,
                              levels = c("astro", "oligo"),   
                              labels = c("IDH-A", "IDH-O")) 
x <- as.character(meta.sub$orig.ident)
col.annot.df$sample <- ifelse(grepl("S[0-9]{4}", x),
                              sub(".*(S[0-9]{4}).*", "\\1", x),
                              NA_character_)

mat.to.plot <- mat.scaled[all.genes.plot, rownames(col.annot.df)]

manual_colors <- c(
  "Inflammatory"         = "#8B1A1A",
  "Microglia"            = "#F4A582",
  "Macrophage"           = "#92C5DE",
  "SPP1+ Macrophage"     = "#08306B",
  "Interferon-response"  = "#BDBDBD"
)

grade_colors <- c("2" = "skyblue", "3" = "royalblue", "4" = "darkblue")
cancer_colors <- c("IDH-A" = "darkgreen", "IDH-O" = "lightgreen")


sample_colors <- c(
  "S0021" = "#1B9E77",
  "S0031" = "#D95F02",
  "S0055" = "#7570B3",
  "S0034" = "#E7298A",
  "S0029" = "#66A61E",
  "S0040" = "#E6AB02",
  "S0028" = "#A6761D",
  "S0001" = "#666666",
  "S0036" = "#1F78B4",
  "S0053" = "#B2DF8A",
  "S0004" = "#FB9A99",
  "S0054" = "#33A02C",
  "S0056" = "#E31A1C",
  "S0035" = "#6A3D9A",
  "S0023" = "#FF7F00",
  "S0018" = "#CAB2D6",
  "S0016" = "#B15928",
  "S0022" = "#8DD3C7",
  "S0020" = "#FFFFB3",
  "S0049" = "#BEBADA",
  "S0060" = "#FB8072",
  "S0037" = "#80B1D3",
  "S0041" = "#FDB462",
  "S0061" = "#B3DE69",
  "S0045" = "#FCCDE5",
  "S0059" = "#BC80BD",
  "S0017" = "#CCEBC5",
  "S0057" = "#FFED6F",
  "S0024" = "#17BECF",
  "S0046" = "#9467BD",
  "S0050" = "#8C564B",
  "S0002" = "#E377C2",
  "S0032" = "#7F7F7F",
  "S0033" = "#BCBD22",
  "S0015" = "#AEC7E8",
  "S0052" = "#FFBB78",
  "S0011" = "#98DF8A",
  "S0048" = "#FF9896",
  "S0030" = "#C5B0D5",
  "S0012" = "#C49C94",
  "S0005" = "#F7B6D2",
  "S0010" = "#DBDB8D",
  "S0013" = "#9EDAE5",
  "S0006" = "#393B79",
  "S0042" = "#637939",
  "S0027" = "#8C6D31",
  "S0044" = "#843C39",
  "S0039" = "#7B4173",
  "S0009" = "#3182BD"
)
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "#F7F7F7", "#B2182B")
)

genes_interet <- c(
  "P2RY12", "CX3CR1", "APBB1IP",      
  "CCL3", "NR4A1", "OLR1",            
  "CD163", "MRC1", "MS4A4A", "F13A1",
  "SPP1","FCGR3A","CPM",
  "IFI44L","MX2","PARP14"
)
genes_interet <- intersect(genes_interet, rownames(mat.to.plot))

row_anno <- rowAnnotation(
  link = anno_mark(
    at = which(rownames(mat.to.plot) %in% genes_interet),
    labels = rownames(mat.to.plot)[rownames(mat.to.plot) %in% genes_interet],
    labels_gp = gpar(fontsize = 16, fontface = "bold"),
    link_width = unit(6, "mm"),
    link_gp = gpar(col = "black", lwd = 0.6)
  )
)

top_anno <- HeatmapAnnotation(
  `Myeloid subset` = col.annot.df$cluster,
  Grade   = col.annot.df$grade,
  Cancer  = col.annot.df$cancer,
  Sample  = col.annot.df$sample,
  
  col = list(
    `Myeloid subset` = manual_colors,
    Grade   = grade_colors,
    Cancer  = cancer_colors,
    Sample  = sample_colors
  ),
  
  show_legend = c(
    `Myeloid subset` = TRUE,
    Grade   = TRUE,
    Cancer  = TRUE,
    Sample  = FALSE
  ),
  
  annotation_name_gp = gpar(fontsize = 15, fontface = "bold"),
  simple_anno_size = unit(0.6, "cm"),
  
  annotation_legend_param = list(
    `Myeloid subset` = list(title = "Myeloid subset",
                            title_gp = gpar(fontsize = 20, fontface = "bold"),
                            labels_gp = gpar(fontsize = 18),
                            grid_height = unit(1, "cm"),
                            grid_width = unit(1, "cm")),
    Grade   = list(title = "Grade",
                   title_gp = gpar(fontsize = 20, fontface = "bold"),
                   labels_gp = gpar(fontsize = 18),
                   grid_height = unit(1, "cm"),
                   grid_width = unit(1, "cm")),
    Cancer  = list(title = "Cancer",
                   title_gp = gpar(fontsize = 20, fontface = "bold"),
                   labels_gp = gpar(fontsize = 18),
                   grid_height = unit(1, "cm"),
                   grid_width = unit(1, "cm"))
  )
)

lgd_expression <- Legend(
  col_fun = col_fun,
  title = "Expression",
  at = c(-2, 0, 2),
  labels = c(-2, "0", 2),
  title_gp = gpar(fontsize = 20, fontface = "bold"),
  labels_gp = gpar(fontsize = 18),
  legend_height = unit(4, "cm"),
  legend_width = unit(1.2, "cm"),
  grid_width = unit(0.7, "cm"),
  title_gap = unit(6, "mm")   
)

ht <- Heatmap(
  mat.to.plot,
  name = "Expression",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  column_split = col.annot.df$cluster,
  cluster_column_slices = FALSE,
  column_gap = unit(1, "mm"),
  column_title = NULL,
  
  show_column_names = FALSE,
  show_row_names = FALSE,
  border = FALSE,
  
  top_annotation = top_anno,
  right_annotation = row_anno,
  
  show_heatmap_legend = FALSE,   
  
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y, w + unit(0.3, "pt"), h + unit(0.3, "pt"),
              gp = gpar(fill = fill, col = NA))
  }
)

## ============================================================
## Fig6.C supp — Heatmap gene expression
## ============================================================
pdf(paste0("output/", "/fig6_C_supp.pdf"), width = 16, height = 10, useDingbats = FALSE)
draw(ht,
     heatmap_legend_list = list(lgd_expression),
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()



################################
# Heatmap gene expression Main #
################################
# rm(list = ls())
# gc()

all.myeloid <- readRDS('all.myeloid.rds')
cells.to.plot <- read.table('heatmap_main_cells_to_plot.tsv', header = FALSE, sep = "\t")[,1]
genes.to.plot <- read.table('heatmap_main_genes_to_plot.tsv', header = FALSE, sep = "\t")[,1]

all.genes.plot <- c(genes.to.plot)

Ncell = 1000
N.per.cat = 20

mat <- all.myeloid@assays$SCT@data
mat.scaled <- t(scale(t(mat)))

all.genes.plot <- c(genes.to.plot)

# Remove cycling cells

cells.to.plot <- cells.to.plot[1:4000]
col.annot.df <- data.frame(
  cluster = c(rep("Microglia", Ncell),
              rep("Inflammatory", Ncell),
              rep("Macrophage", Ncell),
              rep("SPP1+ Macrophage", Ncell))
)
rownames(col.annot.df) <- cells.to.plot

col.annot.df$cluster <- factor(col.annot.df$cluster,
                               levels = c("Microglia", "Inflammatory", "Macrophage","SPP1+ Macrophage"))

all.genes.plot <- intersect(genes.to.plot, rownames(mat.scaled))
cells.to.plot <- intersect(cells.to.plot, colnames(mat.scaled))

mat.to.plot <- mat.scaled[
  all.genes.plot,
  cells.to.plot,
  drop = FALSE
]

col.annot.df <- col.annot.df[cells.to.plot, , drop = FALSE]

manual_colors <- c(
  "Inflammatory" = "#8B1A1A",
  "Microglia"    = "#F4A582",
  "Macrophage"   = "#92C5DE",
  "SPP1+ Macrophage"   = "#08306B"
)

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "#F7F7F7", "#B2182B")
)

# Genes to annotate on the right
genes_interet <- c(
  "P2RY12", "CX3CR1", "APBB1IP",      
  "CCL3", "NR4A1", "OLR1",            
  "CD163", "MRC1", "MS4A4A", "F13A1",
  "SPP1","FCGR3A","CPM"
)
genes_interet <- intersect(genes_interet, rownames(mat.to.plot))

row_anno <- rowAnnotation(
  link = anno_mark(
    at = which(rownames(mat.to.plot) %in% genes_interet),
    labels = rownames(mat.to.plot)[rownames(mat.to.plot) %in% genes_interet],
    labels_gp = gpar(fontsize = 16, fontface = "bold"),
    link_width = unit(6, "mm"),
    link_gp = gpar(col = "black", lwd = 0.6)
  )
)

top_anno <- HeatmapAnnotation(
  'Myeloid subset' = col.annot.df$cluster,
  col = list('Myeloid subset' = manual_colors),
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.4, "cm"),
  annotation_legend_param = list(
    'Myeloid subset' = list(title = "Myeloid subset",
                            title_gp = gpar(fontsize = 20, fontface = "bold"),
                            labels_gp = gpar(fontsize = 18),
                            grid_height = unit(1, "cm"),
                            grid_width = unit(1, "cm"))
  )
)


lgd_expression <- Legend(
  col_fun = col_fun,
  title = "Expression",
  at = c(-2, 0, 2),
  labels = c(-2, "0", 2),
  title_gp = gpar(fontsize = 20, fontface = "bold"),
  labels_gp = gpar(fontsize = 18),
  legend_height = unit(4, "cm"),
  legend_width = unit(1.2, "cm"),
  grid_width = unit(0.7, "cm"),
  title_gap = unit(6, "mm")   
)


ht <- Heatmap(
  mat.to.plot,
  name = "Expression",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  column_split = col.annot.df$cluster,
  cluster_column_slices = FALSE,
  column_gap = unit(1, "mm"),
  column_title = NULL,
  
  show_column_names = FALSE,
  show_row_names = FALSE,
  border = FALSE,
  
  top_annotation = top_anno,
  right_annotation = row_anno,
  
  show_heatmap_legend = FALSE,   
  
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y, w + unit(0.3, "pt"), h + unit(0.3, "pt"),
              gp = gpar(fill = fill, col = NA))
  }
)

## ============================================================
## Fig6.B main — Heatmap gene expression
## ============================================================
pdf(paste0("output/", "/Figure_6_B.pdf"), width = 12, height = 9, useDingbats = FALSE)
draw(ht,
     heatmap_legend_list = list(lgd_expression),
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()


######################
### Proportion plots #
######################

all.myeloid <- readRDS('all.myeloid.filt.rds')
data_all <-  read.table("cell_annotations_minus_S0009.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

ordered_clusters <- c("Microglia", "Inflammatory", "Macrophage", "SPP1+ Macrophage")

custom_colors <- c(
  "Microglia"        = "#F4A582",
  "Inflammatory"     = "#8B1A1A",
  "Macrophage"       = "#92C5DE",
  "SPP1+ Macrophage" = "#08306B"
)

data_all$cluster <- factor(data_all$cluster, levels = ordered_clusters)

meta <- all.myeloid@meta.data
meta$cell <- rownames(meta)
meta$disease <- trimws(as.character(meta$disease))

meta <- meta %>%
  left_join(data_all %>% select(cell.name, cluster), by = c("cell" = "cell.name"))

## ============================================================
## Fig6.F main — Proportion per cells, per grade, per cancer type
## ============================================================
compute_props <- function(meta, disease_filter, grades_keep) {
  meta %>%
    filter(disease == disease_filter, grade %in% grades_keep, !is.na(cluster)) %>%
    mutate(
      grade = factor(grade, levels = grades_keep),
      cluster = factor(cluster, levels = ordered_clusters)
    ) %>%
    count(grade, cluster) %>%
    group_by(grade) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
}

plot_proportions <- function(props_df, title) {
  ggplot(props_df, aes(x = grade, y = prop, fill = cluster)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_fill_manual(values = custom_colors, name = "Myeloid subset", limits = ordered_clusters) +
    scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
    labs(title = title, x = "Grade", y = "Proportion per cell") +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      legend.key.size = unit(1, "cm")
    )
}

props_astro <- compute_props(meta, disease_filter = "astro", grades_keep = c(2, 3, 4))
props_oligo <- compute_props(meta, disease_filter = "oligo", grades_keep = c(2, 3))

p_astro <- plot_proportions(props_astro, "IDH-A")
p_oligo <- plot_proportions(props_oligo, "IDH-O")

p1_stacked <- p_astro + p_oligo + plot_layout(guides = "collect")
ggsave(
  filename = "output/Figure_6_F.pdf",
  plot = p1_stacked,
  width = 8,
  height = 6,
  units = "in"
)

## ============================================================
## Fig6.C main — Proportion per cells, per cancer type
## ============================================================
compute_props_disease <- function(meta) {
  meta %>%
    filter(!is.na(cluster), disease %in% c("astro", "oligo")) %>%
    mutate(
      disease = factor(disease, levels = c("astro", "oligo"),
                       labels = c("IDH-A", "IDH-O")),
      cluster = factor(cluster, levels = ordered_clusters)
    ) %>%
    count(disease, cluster) %>%
    group_by(disease) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
}

props_disease_total <- compute_props_disease(meta)

p_disease_total <- ggplot(props_disease_total, aes(x = disease, y = prop, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) +
  scale_fill_manual(values = custom_colors, name = "Myeloid subset", limits = ordered_clusters) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
  labs(title = "IDH-A vs IDH-O", x = NULL, y = "Proportion per cell") +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm")
  )


ggsave(
  filename = "output/Fig_6_C.pdf",
  plot = p_disease_total,
  width = 8,
  height = 6,
  units = "in"
)


## ============================================================
## Fig6.E main — Proportion per samples, per grade, per cancer type
## ============================================================
compute_props_sample <- function(meta, disease_filter, grades_keep, sample_col = "sample") {
  meta %>%
    filter(disease == disease_filter, grade %in% grades_keep, !is.na(cluster)) %>%
    mutate(
      grade = factor(grade, levels = grades_keep),
      cluster = factor(cluster, levels = ordered_clusters)
    ) %>%
    count(.data[[sample_col]], grade, cluster, .drop = FALSE) %>%
    group_by(.data[[sample_col]], grade) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    filter(!is.nan(prop))
}

get_annotations <- function(props_df, comparisons_list) {
  stat_res <- suppressWarnings(
    compare_means(prop ~ grade, data = props_df, group.by = "cluster", method = "wilcox.test")
  )
  
  keep_pair <- function(g1, g2) {
    any(vapply(comparisons_list, function(cp) setequal(cp, c(g1, g2)), logical(1)))
  }
  stat_res <- stat_res %>% rowwise() %>% filter(keep_pair(group1, group2)) %>% ungroup()
  
  y_max <- props_df %>% group_by(cluster) %>% summarise(ymax = max(prop, na.rm = TRUE), .groups = "drop")
  grade_levels <- levels(props_df$grade)
  
  stat_res %>%
    left_join(y_max, by = "cluster") %>%
    mutate(xmin_pos = match(group1, grade_levels), xmax_pos = match(group2, grade_levels)) %>%
    group_by(cluster) %>%
    arrange(match(paste(group1, group2), sapply(comparisons_list, paste, collapse = " "))) %>%
    mutate(step = row_number(), y_position = ymax * (1.10 + 0.15 * (step - 1))) %>%
    ungroup() %>%
    mutate(label = paste0("p = ", format.pval(p, digits = 2, eps = 0.001)))
}

plot_proportions_sample <- function(props_df, title, comparisons_list, wrap_width = 12) {
  ann <- get_annotations(props_df, comparisons_list)
  
  ggplot(props_df, aes(x = grade, y = prop, fill = cluster)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
    geom_jitter(color = "black", width = 0.15, size = 1.8, alpha = 0.8) +
    geom_signif(
      data = ann,
      aes(xmin = xmin_pos, xmax = xmax_pos, y_position = y_position, annotations = label),
      manual = TRUE, inherit.aes = FALSE, tip_length = 0.01, textsize = 4
    ) +
    facet_wrap(~ cluster, scales = "free_y", nrow = 1,
               labeller = labeller(cluster = label_wrap_gen(width = wrap_width))) +
    scale_fill_manual(values = custom_colors, name = "Myeloid subset", limits = ordered_clusters) +
    scale_y_continuous(labels = function(x) paste0(x * 100, "%"),
                       expand = expansion(mult = c(0.05, 0.35))) +
    labs(title = title, x = "Grade", y = "Proportion per sample") +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 13),
      legend.position = "none",
      strip.text = element_text(size = 15, face = "bold", lineheight = 0.9)
    )
}

props_astro_sample <- compute_props_sample(meta, "astro", c(2, 3, 4))
props_oligo_sample <- compute_props_sample(meta, "oligo", c(2, 3))

comparisons_astro <- list(c("2", "3"), c("2", "4"), c("3", "4"))
comparisons_oligo <- list(c("2", "3"))

p_astro_sample <- plot_proportions_sample(props_astro_sample, "IDH-A", comparisons_astro)
p_oligo_sample <- plot_proportions_sample(props_oligo_sample, "IDH-O", comparisons_oligo)

p2 = p_astro_sample / p_oligo_sample

ggsave(
  filename = "output/Fig_6_E.pdf",
  plot = p2,
  width = 10,
  height = 9,
  units = "in"
)

## ============================================================
## Fig6.F supp — Pair plot ini/rec
## ============================================================

ordered_clusters_paired <- setdiff(ordered_clusters, "SPP1+ Macrophage")

compute_props_patient <- function(meta, disease_filter, patient_col = "patient", status_col = "status") {
  meta %>%
    filter(disease == disease_filter, !is.na(cluster)) %>%
    mutate(cluster = factor(cluster, levels = ordered_clusters), status = .data[[status_col]]) %>%
    count(.data[[patient_col]], status, cluster, .drop = FALSE) %>%
    rename(patient = 1) %>%
    group_by(patient, status) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    filter(!is.nan(prop))
}

keep_paired_patients <- function(props_df) {
  paired_ids <- props_df %>%
    distinct(patient, status) %>%
    count(patient) %>%
    filter(n == 2) %>%
    pull(patient)
  props_df %>% filter(patient %in% paired_ids)
}

get_paired_pvalues <- function(props_df, status_levels = c("ini", "rec1")) {
  wide <- props_df %>%
    select(patient, status, cluster, prop) %>%
    pivot_wider(names_from = status, values_from = prop)
  
  results <- wide %>%
    group_by(cluster) %>%
    summarise(
      p = tryCatch(
        wilcox.test(.data[[status_levels[1]]], .data[[status_levels[2]]], paired = TRUE)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    )
  
  y_max <- props_df %>% group_by(cluster) %>% summarise(ymax = max(prop, na.rm = TRUE), .groups = "drop")
  
  results %>%
    left_join(y_max, by = "cluster") %>%
    mutate(
      label = ifelse(is.na(p), "p = NA", paste0("p = ", format.pval(p, digits = 2, eps = 0.001))),
      y_position = ymax * 1.15
    )
}

plot_paired <- function(props_df, title, status_levels = c("ini", "rec1"), wrap_width = 12) {
  props_df <- props_df %>% mutate(status = factor(status, levels = status_levels))
  ann <- get_paired_pvalues(props_df, status_levels)
  
  ggplot(props_df, aes(x = status, y = prop)) +
    geom_boxplot(aes(fill = cluster), outlier.shape = NA, alpha = 0.6, width = 0.5) +
    geom_line(aes(group = patient), color = "grey60", linewidth = 0.4, alpha = 0.7) +
    geom_point(color = "black", size = 1.8, alpha = 0.8) +
    geom_text(data = ann, aes(x = 1.5, y = y_position, label = label), inherit.aes = FALSE, size = 3) +
    facet_wrap(~ cluster, scales = "free_y", nrow = 1,
               labeller = labeller(cluster = label_wrap_gen(width = wrap_width))) +
    scale_fill_manual(values = custom_colors, name = "Cluster", limits = ordered_clusters) +
    scale_y_continuous(labels = function(x) paste0(x * 100, "%"),
                       expand = expansion(mult = c(0.05, 0.25))) +
    labs(title = title, x = NULL, y = "Proportion per sample") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 13),
      legend.position = "none",
      strip.text = element_text(size = 15, face = "bold", lineheight = 0.9)
    )
}

props_astro_paired <- compute_props_patient(meta, "astro") %>% keep_paired_patients()
props_oligo_paired <- compute_props_patient(meta, "oligo") %>% keep_paired_patients()

p_astro_paired <- plot_paired(props_astro_paired, "IDH-A")
p_oligo_paired <- plot_paired(props_oligo_paired, "IDH-O")

p4 <- p_astro_paired / p_oligo_paired

ggsave(
  filename = "Fig_6_F_supp.pdf",
  plot = p4,
  width = 12,
  height = 9,
  units = "in"
)



## ============================================================
## Fig6.B supp and Fig6.A main — UMAP
## ============================================================
all.myeloid <- readRDS('all.myeloid.filt.rds')
data_all <-  read.table("cell_annotations_minus_S0009.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

meta <- all.myeloid@meta.data
meta$cell <- rownames(meta)

meta <- meta %>%
  left_join(data_all %>% select(cell.name, cluster, umaphar_1, umaphar_2), by = c("cell" = "cell.name"))

meta <- meta %>% filter(!is.na(cluster), !is.na(umaphar_1), !is.na(umaphar_2))

all.myeloid_slim <- subset(all.myeloid, cells = meta$cell)

embeddings <- as.matrix(meta[, c("umaphar_1", "umaphar_2")])
rownames(embeddings) <- meta$cell
colnames(embeddings) <- c("UMAPhar_1", "UMAPhar_2")
embeddings <- embeddings[colnames(all.myeloid_slim), ]

all.myeloid_slim[["UMAP_ICA_state"]] <- CreateDimReducObject(
  embeddings = embeddings,
  key = "UMAPhar_",
  assay = DefaultAssay(all.myeloid_slim)
)

all.myeloid_slim$cluster       <- meta$cluster[match(colnames(all.myeloid_slim), meta$cell)]
all.myeloid_slim$disease       <- all.myeloid_slim@meta.data$disease
all.myeloid_slim$label_program <- all.myeloid_slim@meta.data$label_program
all.myeloid_slim$grade         <- as.character(all.myeloid_slim@meta.data$grade)

# Cancer : valeurs d'origine astro / oligo, affichées IDH-A / IDH-O
all.myeloid_slim$cancer <- factor(trimws(as.character(all.myeloid_slim$disease)),
                                  levels = c("astro", "oligo"),
                                  labels = c("IDH-A", "IDH-O"))
table(all.myeloid_slim$cancer, useNA = "ifany")   # vérification

grade_colors  <- c("2" = "skyblue", "3" = "royalblue", "4" = "darkblue")
cancer_colors <- c("IDH-A" = "darkgreen", "IDH-O" = "lightgreen")

ordered_clusters <- c("Microglia", "Inflammatory", "Macrophage", "SPP1+ Macrophage")
custom_colors <- c(
  "Microglia"        = "#F4A582",
  "Inflammatory"     = "#8B1A1A",
  "Macrophage"       = "#92C5DE",
  "SPP1+ Macrophage" = "#08306B"
)

patient_colors <- c(
  P0014 = "#F8766D", P0016 = "#EF7F49", P0018 = "#E58700", P0025 = "#D89000",
  P0027 = "#C99800", P0029 = "#B79F00", P0002 = "#A3A500", P0004 = "#8AAB00",
  P0006 = "#6BB100", P0008 = "#39B600", P0010 = "#00BA38", P0012 = "#00BD5F",
  P0021 = "#00BF7D", P0023 = "#00C097", P0030 = "#00C0AF", P0015 = "#00BFC4",
  P0017 = "#00BCD8", P0019 = "#00B7E9", P0026 = "#00B0F6", P0028 = "#00A7FF",
  P0001 = "#619CFF", P0003 = "#9590FF", P0005 = "#B983FF", P0007 = "#D376FF",
  P0009 = "#E76BF3", P0011 = "#F564E3", P0020 = "#FD61D1", P0022 = "#FF62BC",
  P0024 = "#FF67A4", P0031 = "#FE6E8A"
)

program_colors <- c(
  "Complement_Immunosuppressive" = "#92C5DE",
  "Scavenger_Immunosuppressive"  = "#08306B",
  "IL1B_CD83_TNF_Inflammatory"   = "#8B1A1A",
  "Inflammatory_microglia"       = "#F4A582",
  "Microglia"                    = "#E69F00",
  "Macrophage"                   = "#1B9E77",
  "Monocyte"                     = "#66A61E",
  "Neutrophil"                   = "#7570B3",
  "cDC"                          = "#E7298A",
  "IFN_Response"                 = "#984EA3",
  "Hypoxia"                      = "#A6761D",
  "HS_UPR"                       = "#17BECF",
  "G1_S"                         = "#999999",
  "G2_M"                         = "#4D4D4D"
)

all.myeloid_slim$cluster <- factor(all.myeloid_slim$cluster, levels = ordered_clusters)

program_levels <- c(
  "Complement_Immunosuppressive",
  "Scavenger_Immunosuppressive",
  "IL1B_CD83_TNF_Inflammatory",
  "Inflammatory_microglia",
  "Microglia",
  "Macrophage",
  "Monocyte",
  "Neutrophil",
  "cDC",
  "IFN_Response",
  "Hypoxia",
  "HS_UPR",
  "G1_S",
  "G2_M"
)

all.myeloid_slim$label_program <- factor(all.myeloid_slim$label_program,
                                         levels = program_levels)

theme_pub <- theme_classic(base_size = 11) +
  theme(
    axis.line        = element_line(linewidth = 0.35, color = "black"),
    axis.text        = element_text(color = "black", size = 10),
    axis.title       = element_text(size = 12, face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 22),
    legend.title     = element_text(size = 18, face = "bold"),
    legend.text      = element_text(size = 12),
    legend.key.size  = unit(0.6, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.background     = element_blank(),
    legend.box.background = element_blank(),
    aspect.ratio     = 1,
    panel.border     = element_blank()
  )

pt_size <- 0.4

make_dimplot <- function(group.by, cols = NULL, title, legend_title = title,
                         shuffle = FALSE, label = FALSE, legend_ncol = 1) {
  DimPlot(all.myeloid_slim, group.by = group.by,
          reduction = "UMAP_ICA_state", pt.size = pt_size,
          cols = cols, shuffle = shuffle, seed = 1,
          label = label, label.size = 3, repel = TRUE) +
    theme_pub +
    ggtitle(title) +
    labs(color = legend_title, x = "UMAP 1", y = "UMAP 2") +
    guides(color = guide_legend(ncol = legend_ncol,
                                override.aes = list(size = 4)))
}

p_cluster  <- make_dimplot("cluster", custom_colors, title = "Myeloid subset")
p_cancer   <- make_dimplot("cancer", cancer_colors, "Disease")
p_program  <- make_dimplot("label_program", program_colors, "Program")
p_grade    <- make_dimplot("grade", grade_colors, "Grade")
p_patient  <- make_dimplot("patient", patient_colors, "Patient", legend_ncol = 2)

p_cluster_shuf <- make_dimplot("cluster", custom_colors, title = "Myeloid subset - shuffled",
                               legend_title = "Myeloid subset", shuffle = TRUE)
p_cancer_shuf  <- make_dimplot("cancer", cancer_colors, "Cancer - shuffled",
                               legend_title = "Cancer", shuffle = TRUE)
p_program_shuf <- make_dimplot("label_program", program_colors, "Program - shuffled",
                               legend_title = "Program", shuffle = TRUE)
p_grade_shuf   <- make_dimplot("grade", grade_colors, "Grade - shuffled",
                               legend_title = "Grade", shuffle = TRUE)
p_patient_shuf <- make_dimplot("patient", patient_colors, "Patient - shuffled",
                               legend_title = "Patient", shuffle = TRUE, legend_ncol = 2)

panel_of <- function(p) p + theme(legend.position = "none")

# Légende seule
leg_of <- function(p) {
  g  <- ggplotGrob(p)
  gb <- g$grobs[grepl("^guide-box", vapply(g$grobs, function(x) x$name, character(1)))]
  gb <- Filter(function(x) !inherits(x, "zeroGrob"), gb)
  wrap_elements(full = gb[[1]])
}

save_plot <- function(p, name, width = 10.5, height = 6) {
  fig <- panel_of(p) + leg_of(p) +
    plot_layout(nrow = 1, widths = unit(c(5.2, 4.4), "in"))
  pdf(file.path("output", paste0("umap_", name, ".pdf")),
      width = width, height = height, useDingbats = FALSE)
  print(fig)
  dev.off()
}

save_plot(p_cluster,  "fig_6_A")
save_plot(p_cancer,   "cancer")
save_plot(p_program,  "miller_et_al_program")
save_plot(p_grade,    "grade")
save_plot(p_patient,  "patient")

save_plot(p_cluster_shuf,  "fig_6_A_shuffled")
save_plot(p_cancer_shuf,   "cancer_shuffled")
save_plot(p_program_shuf,  "program_shuffled")
save_plot(p_grade_shuf,    "grade_shuffled")
save_plot(p_patient_shuf,  "patient_shuffled")

combine3 <- function(pa, pb, pc) {
  panel_of(pa) + leg_of(pa) +
    panel_of(pb) + leg_of(pb) +
    panel_of(pc) + leg_of(pc) +
    plot_layout(nrow = 1,
                widths = unit(c(5.2, 1.5, 5.2, 1.5, 5.2, 4.2), "in"))
}

save_combined <- function(fig, name, width = 25, height = 6) {
  pdf(file.path("output", paste0("umap_", name, ".pdf")),
      width = width, height = height, useDingbats = FALSE)
  print(fig)
  dev.off()
}

fig_side      <- combine3(p_grade, p_cancer, p_program)
fig_side_shuf <- combine3(p_grade_shuf, p_cancer_shuf, p_program_shuf)

save_combined(fig_side,      "fig_6_B_supp")
save_combined(fig_side_shuf, "fig_6_B_supp_shuffled")

###################
### TF analyses ###
###################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})


# CONFIG

seurat_rds <- "all.myeloid.filt.rds"

# AUC scores from SCENIC+
auc_files <- c(
  Astro = "direct_gene_based_AUC_astro.tsv.gz",
  Oligo = "direct_gene_based_AUC_oligo.tsv.gz"
)

merged_name <- "Astro_Oligo_common"

cluster_col <- "cluster"
sample_col  <- "sample"
patient_col <- "patient"


extra_covariates <- c("dataset")
meta_covariates  <- setdiff(extra_covariates, "dataset")

cluster_order <- c("Microglia", "Inflammatory", "Macrophage", "SPP1+ Macrophage")

out_dir <- "TF_heatmaps"

# Exclusion list
exclude_ieg <- c()
# Neural / oligodendrocytaire
exclude_neural <- c("TBR1", "POU6F2", "MEIS2", "NR2F1", "NR2F2", "THRB", "SATB2", "NRL",
                    "FOXP2", "OLIG1", "OLIG2", "ASCL1", "MYRF", "POU3F3", "BNC2",
                    "PKNOX2", "ARNT2", "PLAGL1", "RFX7","ARX", "DLX1", "EMX2", "CUX2",
                    "RFX4","NF1","TCF4","TCF7L2","RFX3","GLI3")

use_ubiquitous_exclusion <- TRUE
exclude_ubiquitous <- c()

exclude_tf <- c(exclude_ieg, exclude_neural,
                if (use_ubiquitous_exclusion) exclude_ubiquitous)
exclude_regex <- "^(ZSCAN|ZBTB|ZKSCAN|ZN|SOX|LHX|NFI)"
exclude_by_ds <- list()
exclude_by_ds[[merged_name]] <- c()

# Selection rules
padj_cutoff <- 0.05
rank_by     <- "effect"   
min_z      <- 0.7   
min_margin <- 0.4
min_effect <- 0     

include_repressors <- FALSE  

n_cores         <- 1      
force_recompute <- FALSE  

# Litterature tf
literature_df <- c("SALL1", "SMAD3", "MEF2C", "TAL1",
                   "NFKB1", "RELA", "NFKB2", "BCL3",
                   "IRF8", "KLF4", "PPARG", "MAF",
                   "SREBF1", "CEBPB","HIF1A",
                   "NR1H3", "NFE2L2")

myeloid_tf <- list(
  lineage_master      = c("SPI1", "CEBPA", "CEBPD", "CEBPE", "RUNX1", "GFI1", "IKZF1",
                          "FLI1", "ETV6", "ELF1", "IRF8", "KLF4", "MAFB", "MAF"),
  macrophage_identity = c("PPARG", "NR1H3", "SPIC", "GATA6", "ZEB2", "ID3", "NFIL3"),
  microglia           = c("SALL1", "MEF2C", "MEF2A", "SMAD3", "TAL1", "MITF", "TFEB",
                          "BHLHE40", "BHLHE41"),
  inflammatory_ifn    = c("NFKB1", "NFKB2", "RELA", "REL", "BCL3", "IRF1", "IRF5", "IRF7",
                          "IRF9", "STAT1", "STAT3", "STAT5",
                          "BATF", "HIF1A", "NFE2L2"),
  TAM_SPP1            = c("SREBF1", "XBP1", "ATF6", "ETS2", "RBPJ", "NFATC2")
)

ieg_flag <- c()

extend_literature <- TRUE
if (extend_literature) {
  literature_df <- unique(c(
    literature_df,
    unlist(myeloid_tf[c("macrophage_identity", "microglia",
                        "inflammatory_ifn", "TAM_SPP1")])
  ))
}


col_fontsize <- 8
gap_mm       <- 0       
lit_gap_mm   <- 4       
chunk_gap_mm <- 3     

# Export 
export_formats <- c("pdf", "png")
export_dpi     <- 600   

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# FONCTIONS

clean_tf_name <- function(x) sub("^([A-Za-z0-9]+).*$", "\\1", x)

pick_variant <- function(tf, raw_names, allow_neg = FALSE){
  hits <- raw_names[clean_tf_name(raw_names) == tf]
  if (!allow_neg) hits <- hits[!grepl("-", hits)]
  if (!length(hits)) return(NA_character_)
  if (length(hits) > 1) {
    ng <- suppressWarnings(as.numeric(sub(".*\\(([0-9]+)g\\).*", "\\1", hits)))
    ng[is.na(ng)] <- 0
    sc <- 2 * grepl("_direct_", hits) + grepl("\\+/\\+", hits)
    hits <- hits[order(-sc, -ng)]
  }
  hits[1]
}

load_auc <- function(path, meta, keep_neg = FALSE){
  if (!file.exists(path)) stop("File not found : ", path)
  tf <- read.table(path, sep = "\t", header = TRUE, row.names = 1,
                   check.names = FALSE)
  if (!keep_neg) tf <- tf[, !grepl("-", colnames(tf)), drop = FALSE]
  auc <- t(as.matrix(tf))                       
  colnames(auc) <- sub("___.*$", "", colnames(auc))
  
  cells <- intersect(colnames(auc), rownames(meta))
  if (!length(cells)) stop("No cells between ", path, " and the Seurat")
  m <- meta[cells, , drop = FALSE]
  keep <- m[[cluster_col]] %in% all_clusters
  cells <- cells[keep]; m <- m[keep, , drop = FALSE]
  
  d <- list(
    auc     = auc[, cells, drop = FALSE],
    cluster = factor(m[[cluster_col]], levels = all_clusters),
    sample  = factor(m[[sample_col]]),
    patient = factor(m[[patient_col]])
  )
  for (cv in meta_covariates) {
    if (!cv %in% colnames(m)) stop("Covariable '", cv, "' not in meta.data")
    d[[cv]] <- m[[cv]]
  }
  tab <- table(d$cluster)
  message(sprintf("  %d same cells | per cluster : %s", length(cells),
                  paste(names(tab), tab, sep = "=", collapse = ", ")))
  if (any(tab < 3)) stop("Cluster with < 3 cells in the dataset : ",
                         paste(names(tab)[tab < 3], collapse = ", "))
  d
}


merge_common_regulons <- function(d1, d2, name1 = "Astro", name2 = "Oligo"){
  raw1 <- rownames(d1$auc); raw2 <- rownames(d2$auc)
  clean1 <- unique(clean_tf_name(raw1))
  clean2 <- unique(clean_tf_name(raw2))
  common_tf <- intersect(clean1, clean2)
  if (!length(common_tf)) stop("Aucun régulon commun entre ", name1, " et ", name2)
  
  best1 <- setNames(vapply(common_tf, pick_variant, character(1), raw_names = raw1), common_tf)
  best2 <- setNames(vapply(common_tf, pick_variant, character(1), raw_names = raw2), common_tf)
  ok <- !is.na(best1) & !is.na(best2)
  common_tf <- common_tf[ok]; best1 <- best1[ok]; best2 <- best2[ok]
  
  a1 <- d1$auc[best1, , drop = FALSE]; rownames(a1) <- common_tf
  a2 <- d2$auc[best2, , drop = FALSE]; rownames(a2) <- common_tf
  
  dup <- intersect(colnames(a1), colnames(a2))
  if (length(dup)) {
    message("  ", length(dup), " Barcodes in the two data -> to stay unique")
    colnames(a1)[colnames(a1) %in% dup] <- paste0(colnames(a1)[colnames(a1) %in% dup], "__", name1)
    colnames(a2)[colnames(a2) %in% dup] <- paste0(colnames(a2)[colnames(a2) %in% dup], "__", name2)
  }
  merged_auc <- cbind(a1, a2)
  
  message(sprintf("  Commun regulon %s/%s : %d keep (on %d and %d) | %d + %d = %d cells",
                  name1, name2, length(common_tf), length(clean1), length(clean2),
                  ncol(a1), ncol(a2), ncol(merged_auc)))
  
  d <- list(
    auc     = merged_auc,
    cluster = factor(c(as.character(d1$cluster), as.character(d2$cluster)), levels = all_clusters),
    sample  = factor(c(as.character(d1$sample),  as.character(d2$sample))),
    patient = factor(c(as.character(d1$patient), as.character(d2$patient))),
    dataset = factor(rep(c(name1, name2), c(ncol(a1), ncol(a2))))
  )
  for (cv in meta_covariates) d[[cv]] <- c(d1[[cv]], d2[[cv]])
  d
}

run_lmer <- function(d, covariates = extra_covariates){
  base <- data.frame(sample = d$sample, patient = d$patient)
  for (cv in covariates) base[[cv]] <- d[[cv]]
  
  cov_term <- if (length(covariates)) paste0(" + ", paste(covariates, collapse = " + ")) else ""
  form <- as.formula(paste0("auc ~ in_cluster", cov_term, " + (1 | patient/sample)"))
  message("  modele : ", deparse(form))
  
  test1 <- function(tf, cl){
    x <- d$auc[tf, ]
    if (sd(x) == 0) return(c(pval = NA_real_, effect = NA_real_))
    df <- base
    df$auc <- x
    df$in_cluster <- factor(ifelse(d$cluster == cl, "in", "out"),
                            levels = c("out", "in"))
    fit <- tryCatch(
      lmer(form, data = df,
           control = lmerControl(check.conv.singular = "ignore")),
      error = function(e) NULL)
    if (is.null(fit)) return(c(pval = NA_real_, effect = NA_real_))
    s <- summary(fit)$coefficients
    if (!"in_clusterin" %in% rownames(s)) return(c(pval = NA_real_, effect = NA_real_))
    c(pval = s["in_clusterin", "Pr(>|t|)"], effect = s["in_clusterin", "Estimate"])
  }
  
  all_tf <- rownames(d$auc)
  out <- lapply(all_clusters, function(cl){
    message("  lmer - cluster : ", cl)
    res <- parallel::mclapply(all_tf, test1, cl = cl, mc.cores = n_cores)
    m <- do.call(rbind, res)
    data.frame(TF = all_tf, cluster = cl, pval = m[, "pval"],
               effect = m[, "effect"], stringsAsFactors = FALSE)
  })
  do.call(rbind, out) %>%
    group_by(cluster) %>%
    mutate(padj = p.adjust(pval, method = "BH")) %>%
    ungroup() %>%
    as.data.frame()
}

analyse_dataset <- function(name, d, covariates = extra_covariates){
  message("== ", name, " ==")
  
  cov_tag <- if (length(covariates)) paste0("_cov-", paste(covariates, collapse = "-")) else ""
  cache <- file.path(out_dir, paste0("lmer_", name, cov_tag, ".rds"))
  spec <- NULL
  if (file.exists(cache) && !force_recompute) {
    spec <- readRDS(cache)
    if (!all(rownames(d$auc) %in% spec$TF) ||
        !all(all_clusters %in% spec$cluster)) {
      message("  cache obsolète -> recalcul"); spec <- NULL
    } else message("  lmer : ", cache)
  }
  if (is.null(spec)) { spec <- run_lmer(d, covariates = covariates); saveRDS(spec, cache) }
  
  mean_auc <- sapply(all_clusters, function(cl)
    rowMeans(d$auc[, d$cluster == cl, drop = FALSE]))
  colnames(mean_auc) <- all_clusters
  z <- t(scale(t(mean_auc)))
  bad <- apply(z, 1, function(x) any(!is.finite(x)))
  if (any(bad)) message("  eRegulons null removed : ", sum(bad))
  
  list(name = name, z = z[!bad, , drop = FALSE], mean_auc = mean_auc,
       spec = spec, n_cells = as.integer(table(d$cluster)[all_clusters]))
}

count_tf <- function(res){
  z <- res$z
  raw <- rownames(z)
  best <- vapply(unique(clean_tf_name(raw)), pick_variant, character(1), raw_names = raw)
  best <- unique(best[!is.na(best)])
  peak <- setNames(colnames(z)[apply(z[best, , drop = FALSE], 1, which.max)], best)
  excl <- c(exclude_tf, exclude_by_ds[[res$name]])
  
  sig <- res$spec %>%
    filter(TF %in% best, !is.na(padj), padj < padj_cutoff,
           effect > 0, effect >= min_effect) %>%
    mutate(tf = clean_tf_name(TF), peak_cl = peak[TF])
  sig$z_cl <- if (nrow(sig)) z[cbind(sig$TF, sig$cluster)] else numeric(0)
  sig$excl <- sig$tf %in% excl |
    (nzchar(exclude_regex) & grepl(exclude_regex, sig$tf))
  
  tab <- do.call(rbind, lapply(all_clusters, function(cl){
    s <- sig[sig$cluster == cl, ]
    p <- s[s$peak_cl == cl, ]
    data.frame(cluster          = cl,
               significatifs    = nrow(s),
               pic_dans_cluster = nrow(p),
               z_ok             = sum(p$z_cl >= min_z),
               z_ok_non_exclus  = sum(p$z_cl >= min_z & !p$excl))
  }))
  
  message(sprintf("\n[%s] TF : %d", res$name, length(best)))
  message(sprintf("[%s] TF significant (padj < %.2f, effect > 0) in >= 1 cluster : %d",
                  res$name, padj_cutoff, length(unique(sig$tf))))
  message(sprintf("[%s] TF candidats (pic in cluster + z >= %.1f) : %d, with %d after exclusion",
                  res$name, min_z, sum(tab$z_ok), sum(tab$z_ok_non_exclus)))
  print(tab, row.names = FALSE)
  write.csv(tab, file.path(out_dir, paste0("TF_counts_", res$name, ".csv")),
            row.names = FALSE)
  
  lab <- function(x) sub("^SPP1\\+ Macrophage$", "Macrophage SPP1+", x)
  by_cl <- function(col) paste0(lab(tab$cluster), ", ", tab[[col]], collapse = "; ")
  txt <- sprintf(
    paste0("Of the %d TFs analyzed, %d were significantly enriched in at least one cluster ",
           "(%d significant TF\u2013cluster associations in total: %s). ",
           "Among these, %d TFs peaked in the tested cluster, and %d also passed the z-score ",
           "threshold (%s). After exclusion of the TFs listed below, %d cluster-enriched TFs ",
           "remained (%s). The %d TFs with the largest effect size per cluster are shown in Fig. [main]."),
    length(best), length(unique(sig$tf)), nrow(sig), by_cl("significatifs"),
    sum(tab$pic_dans_cluster), sum(tab$z_ok), by_cl("z_ok"),
    sum(tab$z_ok_non_exclus), by_cl("z_ok_non_exclus"), 8)
  cat("\n---- TEXTE RESULTATS ----\n", strwrap(txt, 100), sep = "\n")
  writeLines(txt, file.path(out_dir, paste0("results_text_", res$name, ".txt")))
  invisible(tab)
}

select_top <- function(res, n, exclude){
  z <- res$z
  raw_names <- rownames(z)
  
  is_excl <- function(tf){
    tf %in% exclude | (nzchar(exclude_regex) & grepl(exclude_regex, tf))
  }
  
  best <- vapply(unique(clean_tf_name(raw_names)), pick_variant, character(1),
                 raw_names = raw_names)
  best <- unique(best[!is.na(best)])
  peak <- setNames(colnames(z)[apply(z[best, , drop = FALSE], 1, which.max)], best)
  
  cand <- res$spec %>%
    filter(TF %in% best, !is.na(padj), padj < padj_cutoff,
           effect > 0, effect >= min_effect) %>%
    mutate(tf = clean_tf_name(TF), peak_cl = peak[TF]) %>%
    filter(cluster == peak_cl)
  
  cand$z_cl <- z[cbind(cand$TF, cand$cluster)]
  cand$margin <- vapply(seq_len(nrow(cand)), function(i){
    zi <- z[cand$TF[i], colnames(z) != cand$cluster[i]]
    cand$z_cl[i] - max(zi)
  }, numeric(1))
  cand$score <- if (rank_by == "z") cand$z_cl else cand$effect
  
  weak <- cand %>% filter(!(z_cl >= min_z), !is_excl(tf))
  cand <- cand %>% filter(z_cl >= min_z)
  
  hit <- cand %>% filter(is_excl(tf))
  if (nrow(hit))
    message(sprintf("  [%s] removed : %s", res$name,
                    paste0(hit$tf, " (", hit$cluster, ")", collapse = ", ")))
  
  sel <- cand %>%
    filter(!is_excl(tf)) %>%
    group_by(cluster) %>%
    arrange(desc(score), padj, .by_group = TRUE) %>%
    slice_head(n = n) %>%
    ungroup() %>%
    transmute(block = cluster, tf = tf, raw = TF, z_cl = z_cl, margin = margin) %>%
    as.data.frame()
  
  cnt <- table(factor(sel$block, levels = all_clusters))
  message(sprintf("  [%s] select (threshold z>=%.1f) : %s", res$name,
                  min_z,
                  paste(names(cnt), cnt, sep = "=", collapse = " | ")))
  message(sprintf("  [%s] significant but removed (z < %.1f) : %d",
                  res$name, min_z, nrow(weak)))
  message(sprintf("  [%s] median of the select ones : %.2f",
                  res$name, if (nrow(sel)) median(sel$margin, na.rm = TRUE) else NA))
  if (!nrow(sel)) stop("No TF ", res$name,
                       " : decrease min_z")
  
  sel <- sel[order(factor(sel$block, levels = all_clusters), -sel$z_cl), ]
  sel$block <- as.character(sel$block)
  sel$margin <- NULL
  sel
}

select_lit <- function(res, lit, exclude, taken){
  z <- res$z
  raw_names <- rownames(z)
  rows <- lapply(lit, function(tf){
    out <- NULL
    act <- pick_variant(tf, raw_names, allow_neg = FALSE)
    if (!is.na(act))
      out <- data.frame(tf = tf, raw = act, label = tf, stringsAsFactors = FALSE)
    if (include_repressors) {
      neg <- raw_names[clean_tf_name(raw_names) == tf & grepl("-", raw_names)]
      for (r in neg) {
        sg <- regmatches(r, regexpr("[+-]/[+-]", r))
        if (!length(sg)) sg <- "-"
        out <- rbind(out, data.frame(tf = tf, raw = r,
                                     label = paste0(tf, " (", sg, ")"),
                                     stringsAsFactors = FALSE))
      }
    }
    out
  })
  df <- do.call(rbind, rows)
  message("  [", res$name, "] literature : ",
          paste(setdiff(lit, df$tf), collapse = ", "))
  if (is.null(df)) return(NULL)
  
  df <- df[df$raw %in% rownames(z) &
             !df$label %in% c(exclude, taken) & !df$tf %in% exclude, , drop = FALSE]
  if (!nrow(df)) return(NULL)
  
  df$peak <- colnames(z)[apply(z[df$raw, , drop = FALSE], 1, which.max)]
  df$z_cl <- z[cbind(df$raw, df$peak)]
  df <- df[order(factor(df$peak, levels = all_clusters), -df$z_cl), ]
  
  data.frame(block = "Literature", tf = df$label, raw = df$raw,
             z_cl = df$z_cl, peak = df$peak, stringsAsFactors = FALSE)
}

build_panel <- function(res, sel){
  stopifnot(!anyDuplicated(sel$tf))
  z <- res$z[sel$raw, , drop = FALSE]
  rownames(z) <- sel$tf
  key <- paste(res$spec$TF, res$spec$cluster, sep = " || ")
  get_mat <- function(what){
    m <- do.call(cbind, lapply(all_clusters, function(cl)
      res$spec[[what]][match(paste(sel$raw, cl, sep = " || "), key)]))
    dimnames(m) <- list(sel$tf, all_clusters)
    m
  }
  list(name = res$name, z = z, sel = sel, n_cells = res$n_cells,
       mean_auc = res$mean_auc[sel$raw, , drop = FALSE],
       effect = get_mat("effect"), padj = get_mat("padj"))
}

export_panel <- function(p, prefix){
  sd_ <- data.frame(
    TF = p$sel$tf, eRegulon = p$sel$raw, block = p$sel$block,
    setNames(as.data.frame(round(p$mean_auc, 5)), paste0("meanAUC_", all_clusters)),
    setNames(as.data.frame(round(p$z, 3)),        paste0("z_", all_clusters)),
    setNames(as.data.frame(signif(p$effect, 3)),  paste0("effect_", all_clusters)),
    setNames(as.data.frame(signif(p$padj, 3)),    paste0("padj_", all_clusters)),
    check.names = FALSE, row.names = NULL)
  if ("peak" %in% colnames(p$sel)) {
    sd_$peak <- p$sel$peak
  }
  write.csv(sd_, file.path(out_dir, paste0(prefix, "_", p$name, "_source_data.csv")),
            row.names = FALSE)
}

block_gaps_mm <- function(blocks){
  lv <- levels(blocks)
  if (length(lv) < 2) return(0)
  ifelse(lv[-1] == "Literature", lit_gap_mm, gap_mm)
}

col_labels <- function(p, n_mode){
  if (n_mode == "labels")
    paste0(all_clusters, " (n=", format(p$n_cells, big.mark = ",", trim = TRUE), ")")
  else all_clusters
}

make_ht <- function(p, lim, cell_h_mm, cell_w_mm, row_fs, n_mode){
  col_fun <- colorRamp2(c(-lim, 0, lim), c("#08519C", "white", "#A50026"))
  blocks  <- factor(p$sel$block, levels = unique(p$sel$block))
  labels  <- paste0(rownames(p$z), ifelse(sub(" .*$", "", rownames(p$z)) %in% ieg_flag,
                                          " \u2020", ""))
  gaps    <- block_gaps_mm(blocks)
  body_h  <- nrow(p$z) * cell_h_mm + sum(gaps)
  
  top_anno <- NULL
  if (n_mode == "top") {
    n_lab <- paste0("n=", format(p$n_cells, big.mark = ",", trim = TRUE))
    top_anno <- HeatmapAnnotation(
      n = anno_text(n_lab, gp = gpar(fontsize = 6.5, col = "grey30"),
                    rot = 0, just = "center", location = 0.5),
      show_annotation_name = FALSE)
  }
  
  Heatmap(
    p$z, name = "Z-score AUC", col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = blocks, row_gap = unit(gaps, "mm"), row_title = NULL,
    row_labels = labels,
    row_names_gp = gpar(fontsize = row_fs, fontface = "bold"),
    top_annotation = top_anno,
    column_labels = col_labels(p, n_mode),
    column_names_side = "bottom", column_names_rot = 90,
    column_names_gp = gpar(fontsize = col_fontsize, fontface = "bold"),
    column_names_centered = FALSE,
    border = TRUE, rect_gp = gpar(col = NA),
    height = unit(body_h, "mm"),
    width  = unit(ncol(p$z) * cell_w_mm, "mm"),
    show_heatmap_legend = FALSE
  )
}

panel_h_in <- function(p, cell_h_mm, n_mode){
  blocks <- factor(p$sel$block, levels = unique(p$sel$block))
  body   <- (nrow(p$z) * cell_h_mm + sum(block_gaps_mm(blocks))) / 25.4
  top    <- if (n_mode == "top") 0.25 else 0
  body + top + max(nchar(col_labels(p, n_mode))) * col_fontsize * 0.55 / 72 + 0.25
}
panel_w_in <- function(p, cell_w_mm, row_fs){
  max(nchar(rownames(p$z))) * row_fs * 0.72 / 72 + 0.15 + ncol(p$z) * cell_w_mm / 25.4
}

render_figure <- function(hts, widths_in, heights_in, lim, legend, file = NULL,
                          device = "pdf"){
  gap_in <- chunk_gap_mm / 25.4
  leg_w  <- if (legend) 0.95 else 0
  W <- sum(widths_in) + gap_in * (length(hts) - 1) + leg_w + 0.1
  H <- max(heights_in) + 0.1
  col_fun <- colorRamp2(c(-lim, 0, lim), c("#08519C", "white", "#A50026"))
  lgd <- Legend(col_fun = col_fun, title = "Row z-score\nof mean AUC",
                at = c(-lim, 0, lim),
                labels = c(sprintf("%.1f", -lim), "0", sprintf("%.1f", lim)),
                direction = "vertical",
                legend_height = unit(3.5, "cm"), grid_width = unit(5, "mm"),
                title_position = "topleft",
                title_gp = gpar(fontsize = 9, fontface = "bold"),
                labels_gp = gpar(fontsize = 9))
  draw_all <- function(){
    grid.newpage()
    x0 <- 0.05
    for (j in seq_along(hts)) {
      pushViewport(viewport(x = unit(x0, "in"), y = unit(1, "npc"),
                            width = unit(widths_in[j], "in"),
                            height = unit(heights_in[j], "in"),
                            just = c("left", "top")))
      draw(hts[[j]], newpage = FALSE)
      popViewport()
      x0 <- x0 + widths_in[j] + gap_in
    }
    if (legend)
      draw(lgd, x = unit(x0 - gap_in + 0.04, "in"), y = unit(1, "npc") - unit(0.03, "in"),
           just = c("left", "top"))
  }
  if (!is.null(file)) {
    switch(device,
           pdf  = cairo_pdf(paste0(file, ".pdf"), width = W, height = H),
           svg  = svg(paste0(file, ".svg"), width = W, height = H),
           png  = png(paste0(file, ".png"), width = W, height = H, units = "in",
                      res = export_dpi, type = "cairo"),
           tiff = tiff(paste0(file, ".tiff"), width = W, height = H, units = "in",
                       res = export_dpi, compression = "lzw", type = "cairo"))
    draw_all()
    dev.off()
  } else draw_all()
  invisible(c(W = W, H = H))
}

make_figure <- function(res_list, n_top, with_lit, prefix,
                        cell_h_mm, cell_w_mm, row_fs, n_mode = "labels"){
  message("\n##### ", prefix, " (top ", n_top, if (with_lit) " + literature", ") #####")
  panels <- lapply(res_list, function(res){
    excl <- c(exclude_tf, exclude_by_ds[[res$name]])
    sel  <- select_top(res, n_top, excl)
    if (with_lit) {
      lit <- select_lit(res, literature_df, excl, taken = sel$tf)
      sel <- dplyr::bind_rows(sel, lit)
      if (!is.null(lit)) {
        cat("\n=== Littérature [", res$name, "] ===\n", sep = "")
        print(lit[, c("tf","peak")], row.names = FALSE)
      }
    }
    build_panel(res, sel)
  })
  for (p in panels) export_panel(p, prefix)
  
  K   <- length(all_clusters)
  lim <- (K - 1) / sqrt(K)   
  
  hts <- lapply(panels, make_ht, lim = lim, cell_h_mm = cell_h_mm,
                cell_w_mm = cell_w_mm, row_fs = row_fs, n_mode = n_mode)
  hin <- vapply(panels, panel_h_in, numeric(1), cell_h_mm = cell_h_mm, n_mode = n_mode)
  win <- vapply(panels, panel_w_in, numeric(1), cell_w_mm = cell_w_mm, row_fs = row_fs)
  
  for (fmt in export_formats) {
    sz <- render_figure(hts, win, hin, lim, legend = TRUE,
                        file = file.path(out_dir, paste0(prefix, "_with_legend")), device = fmt)
    render_figure(hts, win, hin, lim, legend = FALSE,
                  file = file.path(out_dir, paste0(prefix, "_no_legend")), device = fmt)
  }
  message(sprintf("  figure lenght : %.1f x %.1f cm (L x H)", sz["W"] * 2.54, sz["H"] * 2.54))
  render_figure(hts, win, hin, lim, legend = TRUE)   # aperçu RStudio
  cat("\u2020 = immediate response genes.\n")
  invisible(panels)
}


# Analyses
data <- readRDS(seurat_rds)
all.myeloid <- subset(
  data,
  cells = colnames(data)[
    !grepl("-oligo_P0004_ini_S0009$", colnames(data))
  ]
)

meta_all <- all.myeloid@meta.data
rm(all.myeloid); invisible(gc())

stopifnot(all(c(cluster_col, sample_col, patient_col) %in% colnames(meta_all)))
all_clusters <- if (is.null(cluster_order)) levels(factor(meta_all[[cluster_col]])) else cluster_order
stopifnot(all(all_clusters %in% meta_all[[cluster_col]]))

message("== Chargement Astro ==")
d_astro <- load_auc(auc_files[["Astro"]], meta_all, keep_neg = include_repressors)
message("== Chargement Oligo ==")
d_oligo <- load_auc(auc_files[["Oligo"]], meta_all, keep_neg = include_repressors)
d_merged <- merge_common_regulons(d_astro, d_oligo, name1 = "Astro", name2 = "Oligo")
rm(d_astro, d_oligo); invisible(gc())

res_list <- list()
res_list[[merged_name]] <- analyse_dataset(merged_name, d_merged, covariates = extra_covariates)

for (ds in names(res_list)) count_tf(res_list[[ds]])

for (ds in names(res_list)) {
  sel <- select_top(res_list[[ds]], n = 20, exclude = c())
  cat("\n=====", ds, "=====\n")
  print(split(sel$tf, sel$block))
}


# ============================================================
# Fig6.D main : TOP 8 PER CLUSTER
# ============================================================
panels_main <- make_figure(res_list, n_top = 8, with_lit = FALSE,
                           prefix = "Fig_6_D",
                           cell_h_mm = 5, cell_w_mm = 5,      
                           row_fs = 10, n_mode = "labels")

# ============================================================
# Fig6.E supp : ALL TF + LITTERATURE 
# ============================================================
panels_supp <- make_figure(res_list, n_top = 30, with_lit = TRUE,
                           prefix = "Fig_6_E_supp",
                           cell_h_mm = 2.6, cell_w_mm = 14,   
                           row_fs = 7, n_mode = "top")







#####################################################
### Hallmark pathway over-representation analysis ###
#####################################################

all.myeloid <- readRDS("all.myeloid.filt.rds") 

OUTDIR <- "results_myeloid_ora"
dir.create(OUTDIR, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "tables"),  showWarnings = FALSE)
dir.create(file.path(OUTDIR, "figures"), showWarnings = FALSE)

## Parameters
SPECIES      <- "Homo sapiens"
PADJ_THRESH  <- 0.05
TOP_N        <- 10
MARKER_LOGFC <- 0.5
MARKER_PADJ  <- 0.05
CLUSTER_ORDER <- c("Microglia", "Inflammatory", "Macrophage", "SPP1+ Macrophage")


# Marker identification
Idents(all.myeloid) <- all.myeloid$cluster
all.myeloid <- PrepSCTFindMarkers(all.myeloid)

markers <- FindAllMarkers(
  all.myeloid,
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = MARKER_LOGFC
)

markers <- markers %>%
  filter(p_val_adj < MARKER_PADJ) %>%
  mutate(cluster = trimws(as.character(cluster)))

print(table(markers$cluster))   

write.csv(
  markers,
  file.path(OUTDIR, "tables", "FindAllMarkers_all_clusters.csv"),
  row.names = FALSE
)


# Hallmark gene sets
hallmark_df <- tryCatch(
  msigdbr(species = SPECIES, collection = "H"),   
  error = function(e) msigdbr(species = SPECIES, category = "H")
)
hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

# Universe: all genes of the object
universe <- rownames(all.myeloid)


# ORA function
run_ora <- function(query_genes, pathways, universe){
  query_genes <- intersect(query_genes, universe)
  N <- length(universe)
  n <- length(query_genes)
  
  res <- lapply(names(pathways), function(pw){
    pw_genes <- intersect(pathways[[pw]], universe)
    K <- length(pw_genes)
    k <- length(intersect(query_genes, pw_genes))
    if(K == 0 | k == 0) return(NULL)
    
    pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    
    data.frame(
      pathway = pw,
      k = k,
      K = K,
      n = n,
      N = N,
      fold_enrichment = (k/n)/(K/N),
      pval = pval,
      gene_ratio = paste0(k, "/", n)
    )
  })
  
  res <- bind_rows(res)
  res$padj <- p.adjust(res$pval, method = "BH")
  arrange(res, padj)
}


# Run ORA
ora_all <- list()

for(cl in unique(markers$cluster)){
  genes_cl <- markers %>% filter(cluster == cl) %>% pull(gene)
  tmp <- run_ora(genes_cl, hallmark_sets, universe)
  if(nrow(tmp) == 0) next
  tmp$cluster <- cl
  ora_all[[cl]] <- tmp
}

ora_combined <- bind_rows(ora_all)

write.csv(
  ora_combined,
  file.path(OUTDIR, "tables", "ORA_Hallmark_all_clusters.csv"),
  row.names = FALSE
)


# Plot function
make_ora_barplot <- function(df, top_n = TOP_N, padj_thresh = PADJ_THRESH,
                             cluster_order = CLUSTER_ORDER){
  
  df2 <- df %>%
    mutate(
      cluster = trimws(cluster),
      pathway_short = pathway %>%
        str_remove("^HALLMARK_") %>%
        str_replace_all("_", " "),
      logP = -log10(pmax(padj, 1e-300)),
      significance = factor(if_else(padj < padj_thresh, "Significant", "Not significant"),
                            levels = c("Significant", "Not significant"))
    ) %>%
    group_by(cluster) %>%
    mutate(n_sig = sum(padj < padj_thresh)) %>%
    slice_min(order_by = padj, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(cluster_label = paste0(cluster, "\n(", n_sig, " significant)"))
  
  # Panel order = cluster_order
  lv <- unique(df2$cluster_label[order(match(df2$cluster, cluster_order))])
  df2$cluster_label <- factor(df2$cluster_label, levels = lv)
  
  # Most significant at the TOP (first level = bottom, so sort by decreasing padj)
  df2 <- df2 %>%
    arrange(cluster_label, desc(padj)) %>%
    mutate(pathway_id = paste(cluster_label, pathway_short, sep = "___"))
  df2$pathway_id <- factor(df2$pathway_id, levels = unique(df2$pathway_id))
  
  ggplot(df2, aes(x = logP, y = pathway_id, fill = significance)) +
    geom_col(width = 0.75) +
    scale_fill_manual(
      values = c("Significant" = "#2166AC", "Not significant" = "grey70"),
      name = NULL, drop = FALSE
    ) +
    geom_vline(xintercept = -log10(padj_thresh), linetype = "dashed",
               linewidth = 0.5, colour = "grey30") +
    facet_wrap(~cluster_label, nrow = 1, scales = "free_y") +
    scale_y_discrete(labels = function(x) str_replace(x, regex("^.*___", dotall = TRUE), "")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = expression(-log[10]("adjusted p-value")),
      y = NULL,
      title = "Hallmark pathway over-representation analysis"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 13, face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      axis.text.y = element_text(size = 10, face = "bold", colour = "black"),
      axis.text.x = element_text(size = 11, colour = "black"),
      axis.title.x = element_text(size = 13, face = "bold"),
      panel.spacing = unit(0.8, "cm"),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      plot.margin = margin(10, 15, 10, 10)
    )
}


# ============================================================
# Fig6.A supp : ORA analyses 
# ============================================================

p_ora <- make_ora_barplot(ora_combined)

ggsave(file.path(OUTDIR, "figures", "01_ORA_Hallmark_barplot.pdf"),
       p_ora, width = 17, height = 5.5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "figures", "01_ORA_Hallmark_barplot.png"),
       p_ora, width = 17, height = 5.5, dpi = 300)





