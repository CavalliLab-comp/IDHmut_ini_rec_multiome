library(dplyr)
library(patchwork)
library(ggplot2)
library(Seurat)
library(Signac)
library(readxl)
library(tidyr)
library(ggpubr)
library(rstatix)
library(textshape)
library(compositions)
library(easyCODA)
library(vegan)
library(ComplexHeatmap)
library(circlize)
library(data.table)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


### to change
setwd("/mnt/beegfs/home/alyne/stagein")
fig.dir = '/mnt/beegfs/home/alyne/stageout/figures_myeloid'
dir.create(fig.dir)

#read in full object
all.myeloid = readRDS("/mnt/beegfs/home/alyne/stagein/astro_oligo_data/astro_oligo_dataset_all.rds")


N.per.cat = 20
Ncell = 1000

##########################
#to make main umap/heatmap figure
##########################

#remove strange S0009 oligo samples which has high interferon response
samps.to.rem = c( "oligo_P0004_ini_S0009")

all.myeloid.filt = subset(all.myeloid, subset = !(orig.ident %in% samps.to.rem))

#integrate using Harmony
all.myeloid.filt <- IntegrateLayers(
  object = all.myeloid.filt, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.har",
  verbose = FALSE, assay = "RNA"
)
#NN graph, umap, clustering
all.myeloid.filt <- FindNeighbors(all.myeloid.filt, reduction = "integrated.har", dims = 1:30, graph.name = "har_nn")
all.myeloid.filt <- RunUMAP(all.myeloid.filt, reduction = "integrated.har", dims = 1:30, reduction.name = "umap.har", return.model=T)
res = 0.3
all.myeloid.filt <- FindClusters(all.myeloid.filt, resolution = res, cluster.name = paste0("har_clusters",res), graph = "har_nn", algorithm = 1)

#plot UMAPs showing clusters, subtype
DimPlot(all.myeloid.filt, group.by=c(paste0("har_clusters",res),"disease","patient","grade","orig.ident"), 
        reduction = "umap.har", shuffle = T)
ggsave(paste0(fig.dir,"/umaps_main_filt_myeloid_cells_clusters.pdf"), width=25, height=8)


#plot UMAPs showing scores for Miller programs
FeaturePlot(all.myeloid.filt, features = c("Microglia","Complement_Immunosuppressive","Neutrophil","Scavenger_Immunosuppressive","G2_M","IFN_Response",
                                           "Inflammatory_microglia","Macrophage","IL1B_CD83_TNF_Inflammatory","cDC","G1_S","Hypoxia","HS_UPR","Monocyte","SPP1"),
            reduction = "umap.har")
ggsave(paste0(fig.dir,"/umaps_main_filt_myeloid_cells_miller_scores.pdf"), width=25, height=25)



#############
##code to project sample S0009 and assign myeloid state
##we don't use it further here, but it means we don't have to exclude this sample in the survival analysis
##most abundant state is SPP1+ mac, but they are not expressing SPP1...
#############

#get S0009
query = subset(all.myeloid, subset = (orig.ident %in% samps.to.rem))

#project to dataset with myeloid cells from all other samples
anchors <- FindTransferAnchors(
  reference = all.myeloid.filt,
  query = query,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:30,
  recompute.residuals = FALSE  # <-- Add this line
)
query <- MapQuery(
  anchorset = anchors,
  reference = all.myeloid.filt,
  query = query,
  refdata = list(cluster = paste0("har_clusters",res)),
  reference.reduction = "pca",
  reduction.model = "umap.har"
)



df.most = all.myeloid.filt@meta.data[,c("orig.ident","status","patient","grade","sample","disease",paste0("har_clusters",res))]
df.one = query@meta.data[,c("orig.ident","status","patient","grade","sample","disease","predicted.cluster")]

colnames(df.most) = c("orig.ident","status","patient","grade","sample","disease","cluster")
colnames(df.one) = c("orig.ident","status","patient","grade","sample","disease","cluster")
df.all = rbind(df.most,df.one)

df.all$cluster_annotation = df.all$cluster
df.all$cluster_annotation = as.numeric(as.character(df.all$cluster_annotation))
df.all$cluster_annotation[df.all$cluster_annotation==0] = "Microglia"
df.all$cluster_annotation[df.all$cluster_annotation==1] = "Inflammatory"
df.all$cluster_annotation[df.all$cluster_annotation==2] = "Macrophage"
df.all$cluster_annotation[df.all$cluster_annotation==3] = "Macrophage SPP1+"

cell.types = c("Microglia","Inflammatory","Macrophage","Macrophage SPP1+")
df.all$cluster_annotation = factor(df.all$cluster_annotation,levels=cell.types)



#extract cells from Spp1+ and mac and check what average expression of SPP1 is
mat = all.myeloid.filt@assays$SCT@data
cells.spp1 = rownames(df.most[which(df.most$cluster==3),])
cells.nospp1 = rownames(df.most[which(df.most$cluster==2),])
mean.spp1 = mean(mat["SPP1",cells.spp1])
mean.nospp1 = mean(mat["SPP1",cells.nospp1])
print(paste0("mean expression of SPP1+ is ", mean.spp1, " in SPP1+, and ", mean.nospp1, " in SPP1- cells"))
#compare to S0009 cluster 3 cells
cells.S0009.spp1 = rownames(df.one[which(df.one$cluster==3),])
cells.S0009.nospp1 = rownames(df.one[which(df.one$cluster==2),])
mat.S0009 = query@assays$SCT@data
mean.spp1 = mean(mat.S0009["SPP1",cells.S0009.spp1])
mean.nospp1 = mean(mat.S0009["SPP1",cells.S0009.nospp1])
print(paste0("in S0009, mean expression of SPP1+ is ", mean.spp1, " in SPP1+, and ", mean.nospp1, " in SPP1- cells"))

#write full list of cells, umap coordinates, and cluster assignments
fwrite(df.all, paste0(fig.dir,"/cell_annotations_all_samples.tsv.gz"), row.names = T)



#######
##code for UMAP plots
#######




######
##code for markers
######



######
## code for ORA
######



######
## code for heatmaps
######



######
##code for ini-rec and grade proportion plots, boxplots and barplots
######


######
## statistical test for enrichment in grade 4 astros
######

#statistical test to show enrichment in grade 4 astrocytomas
df = all.myeloid.filt@meta.data[,c("orig.ident","disease","grade","cluster_annotation")]
colnames(df) = c("sample","subtype","grade","type")

#compute proportions of SPP1+ in each sample
df.prop <- df %>% group_by(sample, subtype, grade) %>%
  summarise(n_total = n(), n_target = sum(type == "Macrophage SPP1+"),
    prop = n_target / n_total, .groups = "drop")

#assign to groups, either grade 4 astro or other
df.prop$group <- with(
  df.prop,
  ifelse(subtype == "astro" & grade == 4, "Astro grade 4", "Other")
)
df.prop$group = factor(df.prop$group, levels = c("Other","Astro grade 4"))


fit <- glmmTMB(cbind(n_target, n_total - n_target) ~ group, family = betabinomial(), data = df.prop)
summary(fit)

#plot per grade and subtype
ggplot(df.prop, aes(x = interaction(subtype, grade), y = prop)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.02) +
  theme_classic() + ylab("Proportion of SPP1+ macrophages") + xlab("Subtype.grade")
ggsave(paste0(fig.dir,"/boxplot_proportion_SPP1_mac_per_grade_subtype.pdf"), width = 6, height = 4)
















sample.meta <- df.all %>%
  distinct(orig.ident, status, patient, grade, sample, disease)

# Cell counts, including zero counts
df.prop <- df.all %>%
  count(orig.ident, cluster_annotation, name = "n") %>%
  tidyr::complete(
    orig.ident,
    cluster_annotation = cell.types,
    fill = list(n = 0)
  ) %>%
  left_join(sample.meta, by = "orig.ident") %>%
  group_by(orig.ident) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup()


# df.prop = df.prop[which(df.prop$sample!="S0009"),]
plot.list = list()
for(i in 1:4){
  df.plot = df.prop[which(df.prop$disease=="astro" & df.prop$cluster_annotation==cell.types[i]),]
  plot.list[[i]] = ggplot(df.plot, aes(x = status, y = prop)) + geom_boxplot(outlier.shape=NA) +
    geom_line(aes(group = patient)) + geom_point(size = 2) + labs(y="Proportion",x="",title=cell.types[i])
}
plot.list2 = list()
for(i in 1:4){
  df.plot = df.prop[which(df.prop$disease=="oligo" & df.prop$cluster_annotation==cell.types[i]),]
  plot.list2[[i]] = ggplot(df.plot, aes(x = status, y = prop)) + geom_boxplot(outlier.shape=NA) +
    geom_line(aes(group = patient)) + geom_point(size = 2) + labs(y="Proportion",x="",title=cell.types[i])
}
wrap_plots(c(plot.list,plot.list2),ncol=4)
ggsave(paste0(fig.dir,"/cell_type_proportions_myeloid_ini_rec_paired_all_samples.pdf"), width = 15, height = 8)




plot.list = list()
for(i in 1:4){
  df.plot = df.prop[which(df.prop$disease=="astro" & df.prop$cluster_annotation==cell.types[i]),]
  df.plot$grade = factor(df.plot$grade, levels = c(2,3,4))
  plot.list[[i]] = ggplot(df.plot, aes(x = grade, y = prop)) + geom_boxplot(outlier.shape=NA) + stat_compare_means(
    comparisons = list(c("2", "3"), c("2", "4"), c("3", "4")),
    method = "wilcox.test",
    label = "p.format"
  ) + 
    geom_point(size = 2) + labs(y="Proportion",x="",title=cell.types[i])
}
plot.list2 = list()
for(i in 1:4){
  df.plot = df.prop[which(df.prop$disease=="oligo" & df.prop$cluster_annotation==cell.types[i]),]
  df.plot$grade = factor(df.plot$grade, levels = c(2,3))
  plot.list2[[i]] = ggplot(df.plot, aes(x = grade, y = prop)) + geom_boxplot(outlier.shape=NA) +stat_compare_means(
    comparisons = list(c("2", "3")),
    method = "wilcox.test",
    label = "p.format"
  ) + 
    geom_point(size = 2) + labs(y="Proportion",x="",title=cell.types[i])
}
wrap_plots(c(plot.list,plot.list2),ncol=4)
ggsave(paste0(fig.dir,"/cell_type_proportions_myeloid_grade_all_samples.pdf"), width = 15, height = 8)



#get markers of each cluster
all.myeloid.filt = PrepSCTFindMarkers(all.myeloid.filt)
# all.myeloid.filt$har_clusters0.3[which(all.myeloid.filt$har_clusters0.3==4)] = 2
markers = FindAllMarkers(all.myeloid.filt, group.by = paste0("har_clusters",res))

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
fwrite(list(cells = cells.to.plot),paste0(fig.dir,"/heatmap_main_cells_to_plot.tsv.gz"))
fwrite(list(cells = genes.to.plot),paste0(fig.dir,"/heatmap_main_genes_to_plot.tsv.gz"))

#####heatmap
mat = all.myeloid@assays$SCT@data
mat.scaled <- t(scale(t(mat)))
all.cells.plot = c(cells.to.plot)
all.genes.plot = c(genes.to.plot)
mat.to.plot = mat.scaled[all.genes.plot,all.cells.plot]
col.annot.df = data.frame(cluster = c(rep("Microglia",Ncell), rep("Inflammatory",Ncell), rep("Macrophage",Ncell), rep("Complement",Ncell)))

# col.annot.df = data.frame(cluster = c(rep("Microglia",Ncell), rep("Inflammatory",Ncell), rep("Macrophage",Ncell), rep("Cycling",length(cells.to.plot.list[[4]]))))
col.annot.df$cluster = factor(col.annot.df$cluster, levels = c("Microglia","Inflammatory","Macrophage","Complement"))
rownames(col.annot.df) = all.cells.plot
cols.cluster = gg_color_hue(length(levels(col.annot.df$cluster)))
names(cols.cluster) = levels(col.annot.df$cluster)
cols.list.col = list(cluster = cols.cluster)
top_anno <- HeatmapAnnotation(df = col.annot.df, col = cols.list.col,
                              annotation_legend_param = list(cluster = list(labels_gp = gpar(fontsize = 20),
                                                                            title_gp  = gpar(fontsize = 22))))

col_fun <- colorRamp2(c(-2, 0, 2), c("darkblue", "white", "darkred"))
ht <- Heatmap(
  mat.to.plot,
  col = col_fun,
  name = "Gene expression",
  cluster_rows = F,
  cluster_columns = F,
  top_annotation  = top_anno,
  show_row_names = T,
  show_column_names = F,
  row_title = NULL,
  column_title = NULL,
  show_column_dend = F,
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect( x, y, w + unit(0.3, "pt"), h + unit(0.3, "pt"),
               gp = gpar(fill = fill, col = NA) )
  },  heatmap_legend_param = list(legend_height = unit(6, "cm"), labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 22))
)
pdf(paste0(fig.dir,"/heatmap_gex_myeloid_main_filt_nocluster.pdf"), width = 20, height = 20, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()




ht <- Heatmap(
  mat.to.plot,
  col = col_fun,
  name = "Gene expression",
  cluster_rows = T,
  cluster_columns = T,
  top_annotation  = top_anno,
  show_row_names = T,
  show_column_names = F,
  row_title = NULL,
  column_title = NULL,
  show_column_dend = F,
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect( x, y, w + unit(0.3, "pt"), h + unit(0.3, "pt"),
               gp = gpar(fill = fill, col = NA) )
  },  heatmap_legend_param = list(legend_height = unit(6, "cm"), labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 22))
)
pdf(paste0(fig.dir,"/heatmap_gex_myeloid_main_filt_clustered.pdf"), width = 20, height = 20, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()



#write full list of cells, umap coordinates, and cluster assignments
all.myeloid.filt$cluster_annotation = all.myeloid.filt$har_clusters0.3
all.myeloid.filt$cluster_annotation = as.numeric(as.character(all.myeloid.filt$cluster_annotation))
all.myeloid.filt$cluster_annotation[all.myeloid.filt$cluster_annotation==0] = "Microglia"
all.myeloid.filt$cluster_annotation[all.myeloid.filt$cluster_annotation==1] = "Inflammatory"
all.myeloid.filt$cluster_annotation[all.myeloid.filt$cluster_annotation==2] = "Macrophage"
all.myeloid.filt$cluster_annotation[all.myeloid.filt$cluster_annotation==3] = "Macrophage SPP1+"
df.to.write = cbind(data.frame(cell.name = colnames(all.myeloid.filt),
                         cluster = all.myeloid.filt$cluster_annotation),
                    Embeddings(all.myeloid.filt, reduction = "umap.har"))
fwrite(df.to.write,paste0(fig.dir,"/cell_annotations_minus_S0009.tsv.gz"))



#oligo_P0004_ini_S0009
#get markers of cells in sample S0009 compared to other oligo myeloid cells
all.myeloid.oligo = subset(all.myeloid, subset = disease=="oligo")
all.myeloid.oligo = PrepSCTFindMarkers(all.myeloid.oligo)
markers.09 = FindMarkers(all.myeloid.oligo, group.by = "orig.ident", ident.1 = "oligo_P0004_ini_S0009")
genes09.to.plot = rownames(head(markers.09[markers.09$avg_log2FC>0,],N.per.cat+7))
#remove genes already in marker set
genes09.to.plot = genes09.to.plot[!(genes09.to.plot %in% c(genes.to.plot))]
cells.09 = names(all.myeloid.oligo$orig.ident[which(all.myeloid.oligo$orig.ident=="oligo_P0004_ini_S0009")])
cells.09.to.plot = cells.09[sample(length(cells.09),min(Ncell,length(cells.09)))]






mat = all.myeloid@assays$SCT@data
mat.scaled <- t(scale(t(mat)))
all.cells.plot = c(cells.to.plot,cells.09.to.plot)
all.genes.plot = c(genes.to.plot,genes09.to.plot)
mat.to.plot = mat.scaled[all.genes.plot,all.cells.plot]

col.annot.df = data.frame(cluster = c(rep("Microglia",Ncell), rep("Inflammatory",Ncell), rep("Macrophage",Ncell), rep("SPP1+ Macrophage",Ncell),
                                      rep("IFN-responsive",length(cells.09.to.plot))))
col.annot.df$cluster = factor(col.annot.df$cluster, levels = c("Microglia","Inflammatory","Macrophage", "SPP1+ Macrophage",
                                                               "IFN-responsive"))
rownames(col.annot.df) = all.cells.plot
cols.cluster = gg_color_hue(length(levels(col.annot.df$cluster)))
names(cols.cluster) = levels(col.annot.df$cluster)
cols.list.col = list(cluster = cols.cluster)
top_anno <- HeatmapAnnotation(df = col.annot.df, col = cols.list.col,
                              annotation_legend_param = list(cluster = list(labels_gp = gpar(fontsize = 20),
                                                                            title_gp  = gpar(fontsize = 22))))

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht <- Heatmap(
  mat.to.plot,
  col = col_fun,
  name = "Gene expression",
  cluster_rows = F,
  cluster_columns = F,
  top_annotation  = top_anno,
  show_row_names = T,
  show_column_names = F,
  row_title = NULL,
  column_title = NULL,
  show_column_dend = F,
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect( x, y, w + unit(0.3, "pt"), h + unit(0.3, "pt"),
               gp = gpar(fill = fill, col = NA) )
  },  heatmap_legend_param = list(legend_height = unit(6, "cm"), labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 22))
)
pdf(paste0(fig.dir,"/heatmap_gex_myeloid_full_nocluster.pdf"), width = 20, height = 20, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()



#do clustered version just to see how it looks
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht <- Heatmap(
  mat.to.plot,
  col = col_fun,
  name = "Gene expression",
  cluster_rows = T,
  cluster_columns = T,
  top_annotation  = top_anno,
  show_row_names = T,
  show_column_names = F,
  row_title = NULL,
  column_title = NULL,
  show_column_dend = F,
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect( x, y, w + unit(0.3, "pt"), h + unit(0.3, "pt"),
               gp = gpar(fill = fill, col = NA) )
  },
  heatmap_legend_param = list(legend_height = unit(6, "cm"), labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 22))
)
pdf(paste0(fig.dir,"/heatmap_gex_myeloid_full_clustered.pdf"), width = 20, height = 20, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


#write lists of cells and genes for figure
fwrite(list(cells = all.cells.plot),paste0(fig.dir,"/heatmap_full_cells_to_plot.tsv.gz"))
fwrite(list(cells = all.genes.plot),paste0(fig.dir,"/heatmap_full_genes_to_plot.tsv.gz"))







#violin plot of Miller interferon expression score by sample
df <- FetchData(all.myeloid, vars = c("orig.ident", "IFN_Response"))

ggplot(df, aes(x = orig.ident, y = IFN_Response)) +
  geom_violin(trim = FALSE) + geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(fig.dir,"/violin_IFN_score_per_sample.pdf"), width = 15, height = 5)




IFN.sig = rownames(head(markers.09[markers.09$avg_log2FC>0,],100))
library(UCell)
all.myeloid <- AddModuleScore_UCell(all.myeloid,features = list(IFN_Macrophage = IFN.sig))

df <- FetchData(all.myeloid, vars = c("orig.ident", "IFN_Macrophage_UCell","grade"))
df$grade = factor(df$grade)
cols.grade = c("2" = "skyblue", "3" = "royalblue", "4" = "darkblue")
ggplot(df, aes(x = orig.ident, y = IFN_Macrophage_UCell, fill = grade)) +
  geom_boxplot(outlier.shape = NA) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(t = 5, r = 5, b = 30, l = 35)) + scale_fill_manual(values = cols.grade)
ggsave(paste0(fig.dir,"/boxplot_my_IFN_score_per_sample.pdf"), width = 10, height = 5)


signature_core <- c("IFI44L","IFI44","IFI6","IFIT1","IFIT2","IFIT3",
  "ISG15","MX1","MX2","RSAD2","OAS2","OAS3",
  "HERC5","HERC6","CMPK2","USP18","STAT1",
  "SIGLEC1","EPSTI1","XAF1","F13A1","MAFB","PARP9",
  "PARP14","SAMD9L")
all.myeloid <- AddModuleScore_UCell(all.myeloid,features = list(IFN_core_Macrophage = signature_core))

df <- FetchData(all.myeloid, vars = c("orig.ident", "IFN_core_Macrophage_UCell","grade"))
df$grade = factor(df$grade)
cols.grade = c("2" = "skyblue", "3" = "royalblue", "4" = "darkblue")
ggplot(df, aes(x = orig.ident, y = IFN_core_Macrophage_UCell, fill = grade)) +
  geom_boxplot(outlier.shape = NA) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(t = 5, r = 5, b = 30, l = 35)) + scale_fill_manual(values = cols.grade)
ggsave(paste0(fig.dir,"/boxplot_my_core_IFN_score_per_sample.pdf"), width = 10, height = 5)










#list of ligands to check
mac.ligs = c("TLR4","VSIR","SPP1", #from interactions never seen in healthy brain
             "ITGA9","ITGB1","PTPRM","ITGB5","ITGAV","WNT5A",
             "PDGFC","PDGFB","HLA-E","SLC1A3",
             "SORL1","CD74","IL18","TLR2")

lapply(mac.ligs, function(x) {markers.filt[markers.filt$gene==x,]})


#macrophage genes for dot plot from astro/opc specific figures
dot.genes.mac = c("SPP1","ITGB1","ITGAV","ITGA9","ITGB5","PTPRM",
                  "WNT5A", #came out as related to grade 4
                  "TLR4","PDGFB","PDGFC","HLA-E","ABCA1","SLC1A3","CD74","ADAM9","ADAM17",
                  "SORL1","VSIR","IL18")


library(scCustomize)

myeloid.astro = subset(all.myeloid.filt, subset = disease == "astro")
myeloid.oligo = subset(all.myeloid.filt, subset = disease == "oligo")

DotPlot(all.myeloid.filt, features = dot.genes.mac, group.by = "cluster_annotation", scale=T) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+ scale_colour_viridis_c(option = "magma", direction=-1)
ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_scaled.pdf"), width = 8, height = 2.7)

DotPlot(all.myeloid.filt, features = dot.genes.mac, group.by = "cluster_annotation", scale=F) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_colour_viridis_c(option = "magma", direction=-1,limits = c(0, 3.5), oob = scales::squish)
ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_not_scaled.pdf"), width = 8, height = 2.7)


Idents(myeloid.astro) <- "cluster_annotation"
Idents(myeloid.oligo) <- "cluster_annotation"
p1 = DotPlot(myeloid.astro, features = dot.genes.mac, group.by = "cluster_annotation",
             col.min = -0.5, col.max = 1.5, dot.scale = 8)+ scale_colour_viridis_c(option = "magma", direction=-1) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p2 = DotPlot(myeloid.oligo, features = dot.genes.mac, group.by = "cluster_annotation",
             col.min = -0.5, col.max = 1.5, dot.scale = 8)+ scale_colour_viridis_c(option = "magma", direction=-1) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p1 + p2
ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_by_cancer_scaled.pdf"), width = 20, height = 3.5)

p1 = DotPlot(myeloid.astro, features = dot.genes.mac, group.by = "cluster_annotation",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale = F) + scale_colour_viridis_c(option = "magma", direction=-1)+
  theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p2 = DotPlot(myeloid.oligo, features = dot.genes.mac, group.by = "cluster_annotation",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F)+ scale_colour_viridis_c(option = "magma", direction=-1) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p1 + p2
ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_by_cancer_not_scaled.pdf"), width = 20, height = 3.5)





p1 = DotPlot(myeloid.astro, features = dot.genes.mac, group.by = "grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100))+ scale_colour_viridis_c(option = "magma", direction=-1) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p2 = DotPlot(myeloid.oligo, features = dot.genes.mac, group.by = "grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100))+ scale_colour_viridis_c(option = "magma", direction=-1) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p1 + p2
ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_by_cancer_by_grade_scaled.pdf"), width = 17, height = 3.5)

p1 = DotPlot(myeloid.astro, features = dot.genes.mac, group.by = "grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100))+ scale_colour_viridis_c(option = "magma", direction=-1) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p2 = DotPlot(myeloid.oligo, features = dot.genes.mac, group.by = "grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100))+ scale_colour_viridis_c(option = "magma", direction=-1) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p1 + p2
ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_by_cancer_by_grade_not_scaled.pdf"), width = 17, height = 3.5)



dot.genes.mac.final = c(c("SPP1","ITGB1","ITGAV","ITGA9","ITGB5","PDGFB","PDGFC","HLA-E","TLR4","ABCA1","CD74"))
#combined dot plot with key LR
DotPlot(all.myeloid.filt, features = dot.genes.mac.final, group.by = "cluster_annotation",
        col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  scale_colour_viridis_c(option = "magma", direction=-1) + labs(y="",x="") + theme_bw() + scale_size(range = c(0, 6), limits = c(20, 90)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text=element_text(size=12))
ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_combined_scaled.pdf"), width = 10, height = 3.8)

#combined dot plot with key LR
DotPlot(all.myeloid.filt, features = dot.genes.mac.final, group.by = "grade",
        col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  scale_colour_viridis_c(option = "magma", direction=-1) + labs(y="Grade",x="") + theme_bw() + scale_size(range = c(0, 6), limits = c(25, 100)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text=element_text(size=12))
ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_combined_by_grade_scaled.pdf"), width = 8, height = 3.8)


#also make separate ones per cancer
p1 = DotPlot(myeloid.astro, features = dot.genes.mac.final, group.by = "cluster_annotation",
        col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(20, 90))+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 12))+#, axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())+
  scale_colour_viridis_c(option = "magma", direction=-1) +guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2)) +
  labs(y = "", x = "", colour = "Average Expression", size = "Percent Expressed")
p2 = DotPlot(myeloid.oligo, features = dot.genes.mac.final, group.by = "cluster_annotation",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(20, 90)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_colour_viridis_c(option = "magma", direction=-1) + theme(axis.text=element_text(size=12))+guides(colour = guide_colourbar(order = 1),size = guide_legend(order = 2)) +
  labs(y = "",x = "",colour = "Average Expression",size = "Percent Expressed")
p1+p2
ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_by_cancer_scaled.pdf"), width = 18, height = 3.8)



Idents(astro) <- "ICA_mes"
Idents(oligo) <- "ICA_mes"
p1 = DotPlot(astro, features = dot.genes.malig, group.by = "ICA_mes",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(20, 90))+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 12))+#, axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())+
  scale_colour_viridis_c(option = "magma", direction=-1) +guides(
    colour = guide_colourbar(order = 1),size = guide_legend(order = 2)) +
  labs(y = "",x = "",colour = "Average Expression",size = "Percent Expressed")
p2 = DotPlot(oligo, features = dot.genes.malig, group.by = "ICA_mes",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(20, 90)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_colour_viridis_c(option = "magma", direction=-1) + theme(axis.text=element_text(size=12))+guides(
    colour = guide_colourbar(order = 1),size = guide_legend(order = 2)) +
  labs(y = "",x = "",colour = "Average Expression",size = "Percent Expressed")
p1 + p2 + plot_layout(ncol = 2)
ggsave(paste0(fig.dir,"/dotplot_all_malig_cells_mes_by_cancer_scaled.pdf"), width = 18, height = 3.8)




#combined dot plot with key LR
p1 = DotPlot(myeloid.astro, features = dot.genes.mac.final, group.by = "grade",
        col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  scale_colour_viridis_c(option = "magma", direction=-1) + labs(y="Grade",x="") + theme_bw() + scale_size(range = c(0, 6), limits = c(25, 100)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text=element_text(size=12))
p2 = DotPlot(myeloid.oligo, features = dot.genes.mac.final, group.by = "grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  scale_colour_viridis_c(option = "magma", direction=-1) + labs(y="Grade",x="") + theme_bw() + scale_size(range = c(0, 6), limits = c(25, 100)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text=element_text(size=12))
p1+p2
ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_by_cancer_by_grade_scaled.pdf"), width = 18, height = 3.8)


#combined dot plot with key LR
DotPlot(all.myeloid.filt, features = dot.genes.mac.final, group.by = "cluster_annotation",
        col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  scale_colour_viridis_c(option = "magma", direction=-1) + labs(y="",x="") + theme_bw() + scale_size(range = c(0, 6), limits = c(25, 100)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text=element_text(size=12))
ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_combined_not_scaled.pdf"), width = 8, height = 3.8)

#combined dot plot with key LR
DotPlot(all.myeloid.filt, features = dot.genes.mac.final, group.by = "grade",
        col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  scale_colour_viridis_c(option = "magma", direction=-1) + labs(y="Grade",x="") + theme_bw() + scale_size(range = c(0, 6), limits = c(25, 100)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text=element_text(size=12))
ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_combined_by_grade_not_scaled.pdf"), width = 8, height = 3.8)




#for soroban
cancer = "oligo"
data.dir = paste0("/mnt/beegfs/home/alyne/stagein/scenicplus_",cancer)
setwd(data.dir)

library(readxl)

if(cancer=="astro"){
  cor.thresh = 0.45
}else{
  cor.thresh = 0.5
}
res = 0.18

fig.dir = paste0("/mnt/beegfs/home/alyne/stageout/figures_sICA_",cancer,"_",cor.thresh)
dir.create(fig.dir)

#ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#set to "malig" if want to run on malignant cells
fig.str = "malig"
oligo = readRDS(paste0(data.dir,"/oligo.ss_integrated_RPCA_",fig.str,"_",cancer,".rds"))


#read in annotated sICA components
df.annot = as.data.frame(read_xlsx(paste0(fig.dir,"/top_gene_med_",cancer,"_annot.xlsx")))
colnames(df.annot) = df.annot[1,]
df.annot = df.annot[-1,]
prog.orig.order.oligo = colnames(df.annot)



cancer = "astro"
data.dir = paste0("/mnt/beegfs/home/alyne/stagein/scenicplus_",cancer)
if(cancer=="astro"){
  cor.thresh = 0.45
}else{
  cor.thresh = 0.5
}
res = 0.18

fig.dir = paste0("/mnt/beegfs/home/alyne/stageout/figures_sICA_",cancer,"_",cor.thresh)
astro = readRDS(paste0(data.dir,"/oligo.ss_integrated_RPCA_",fig.str,"_",cancer,"_doub_rem.rds"))

#read in annotated sICA components
df.annot = as.data.frame(read_xlsx(paste0(fig.dir,"/top_gene_med_",cancer,"_annot.xlsx")))
colnames(df.annot) = df.annot[1,]
df.annot = df.annot[-1,]
prog.orig.order.astro = colnames(df.annot)

#find top scores for each cell
oligo$top_ICA_oligo = sapply(seq(length(colnames(oligo))), function(x) {prog.orig.order.oligo[which.max(oligo@assays$ICA@layers$data[,x])]})
astro$top_ICA_astro = sapply(seq(length(colnames(astro))), function(x) {prog.orig.order.astro[which.max(astro@assays$ICA@layers$data[,x])]})

oligo$ICA_mes = oligo$ICA_score_states
oligo$ICA_mes[which(oligo$top_ICA_oligo=="Reactive astrocyte")] = "reactive"

astro$ICA_mes = astro$ICA_score_states
astro$ICA_mes[which(astro$top_ICA_astro=="Reactive astrocyte")] = "reactive"

oligo$ICA_mes_grade = paste("grade",oligo$grade,oligo$ICA_mes)
astro$ICA_mes_grade = paste("grade",astro$grade,astro$ICA_mes)

astro$ICA_mes_grade = factor(astro$ICA_mes_grade, levels = c(c("grade 2 astro-like","grade 3 astro-like","grade 4 astro-like",
                                                               "grade 2 reactive","grade 3 reactive","grade 4 reactive"),
                                                             unique(astro$ICA_mes_grade)[!(unique(astro$ICA_mes_grade) %in% c("grade 2 astro-like","grade 3 astro-like","grade 4 astro-like",
                                                                                                                              "grade 2 reactive","grade 3 reactive","grade 4 reactive"))]))
oligo$ICA_mes_grade = factor(oligo$ICA_mes_grade, levels = c(c("grade 2 astro-like","grade 3 astro-like",
                                                               "grade 2 reactive","grade 3 reactive"),
                                                             unique(oligo$ICA_mes_grade)[!(unique(oligo$ICA_mes_grade) %in% c("grade 2 astro-like","grade 3 astro-like",
                                                                                                                              "grade 2 reactive","grade 3 reactive"))]))

dot.genes.astro = c("CD44","TNC","FN1","LAMA4","LTBP1","PTPRM","FZD3","FZD5","LRP5","LRP6")


fig.dir = '/mnt/beegfs/home/alyne/stageout/figures_myeloid'
Idents(astro) <- "ICA_mes_grade"
Idents(oligo) <- "ICA_mes_grade"
p1 = DotPlot(astro, features = dot.genes.astro, group.by = "ICA_mes_grade", idents = c("grade 2 astro-like","grade 3 astro-like","grade 4 astro-like",
                                                                                       "grade 2 reactive","grade 3 reactive","grade 4 reactive"),
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p2 = DotPlot(oligo, features = dot.genes.astro, group.by = "ICA_mes_grade", idents = c("grade 2 astro-like","grade 3 astro-like",
                                                                                       "grade 2 reactive","grade 3 reactive"),
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p1 + p2 #+ plot_layout(guides = "collect")
ggsave(paste0(fig.dir,"/dotplot_astrolike_malig_cells_by_cancer_by_grade_scaled.pdf"), width = 15, height = 3.5)

p1 = DotPlot(astro, features = dot.genes.astro, group.by = "ICA_mes_grade", idents = c("grade 2 astro-like","grade 3 astro-like","grade 4 astro-like",
                                                                                       "grade 2 reactive","grade 3 reactive","grade 4 reactive"),
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p2 = DotPlot(oligo, features = dot.genes.astro, group.by = "ICA_mes_grade", idents = c("grade 2 astro-like","grade 3 astro-like",
                                                                                       "grade 2 reactive","grade 3 reactive"),
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p1 + p2 #+ plot_layout(guides = "collect")
ggsave(paste0(fig.dir,"/dotplot_astrolike_malig_cells_by_cancer_by_grade_not_scaled.pdf"), width = 15, height = 3.5)



Idents(astro) <- "ICA_mes"
Idents(oligo) <- "ICA_mes"

p1 = DotPlot(astro, features = dot.genes.astro, group.by = "ICA_mes", idents = c("reactive","astro-like"),
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p2 = DotPlot(oligo, features = dot.genes.astro, group.by = "ICA_mes", idents = c("reactive","astro-like"),
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p1 + p2 #+ plot_layout(guides = "collect")
ggsave(paste0(fig.dir,"/dotplot_astrolike_malig_cells_by_cancer_scaled.pdf"), width = 15, height = 3.5)

p1 = DotPlot(astro, features = dot.genes.astro, group.by = "ICA_mes", idents = c("reactive","astro-like"),
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p2 = DotPlot(oligo, features = dot.genes.astro, group.by = "ICA_mes", idents = c("reactive","astro-like"),
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p1 + p2 #+ plot_layout(guides = "collect")
ggsave(paste0(fig.dir,"/dotplot_astrolike_malig_cells_by_cancer_not_scaled.pdf"), width = 15, height = 3.5)




dot.genes.malig = c("CD44","TNC","FN1","LAMA4","LTBP1","PDGFRA","KLRC2","ANGPTL2","RIMS2","MIF")
p1 = DotPlot(astro, features = dot.genes.malig, group.by = "ICA_score_states",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p2 = DotPlot(oligo, features = dot.genes.malig, group.by = "ICA_score_states",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p1 + p2 #+ plot_layout(guides = "collect")
ggsave(paste0(fig.dir,"/dotplot_all_malig_cells_by_cancer_by_ICA_score_state_scaled.pdf"), width = 18, height = 3.5)

p1 = DotPlot(astro, features = dot.genes.malig, group.by = "ICA_score_states",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p2 = DotPlot(oligo, features = dot.genes.malig, group.by = "ICA_score_states",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p1 + p2 #+ plot_layout(guides = "collect")
ggsave(paste0(fig.dir,"/dotplot_all_malig_cells_by_cancer_by_ICA_score_state_not_scaled.pdf"), width = 18, height = 3.5)




astro$ICA_mes = factor(astro$ICA_mes, levels = c("astro-like","reactive","OPC-like","neuro-like","intermediate"))
oligo$ICA_mes = factor(oligo$ICA_mes, levels = c("astro-like","reactive","OPC-like","intermediate"))

Idents(astro) <- "ICA_mes"
Idents(oligo) <- "ICA_mes"
p1 = DotPlot(astro, features = dot.genes.malig, group.by = "ICA_mes",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(25, 100))+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 12))+#, axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())+
  scale_colour_viridis_c(option = "magma", direction=-1) +guides(
    colour = guide_colourbar(order = 1),
    size = guide_legend(order = 2)
  ) +
  labs(
    y = "",
    x = "",
    colour = "Average Expression",
    size = "Percent Expressed"
  )
p2 = DotPlot(oligo, features = dot.genes.malig, group.by = "ICA_mes",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(25, 100)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_colour_viridis_c(option = "magma", direction=-1) + theme(axis.text=element_text(size=12))+guides(
    colour = guide_colourbar(order = 1),
    size = guide_legend(order = 2)
  ) +
  labs(
    y = "",
    x = "",
    colour = "Average Expression",
    size = "Percent Expressed"
  )
p1 + p2 + plot_layout(ncol = 2)
ggsave(paste0(fig.dir,"/dotplot_all_malig_cells_mes_by_cancer_scaled.pdf"), width = 18, height = 3.8)


p1 = DotPlot(astro, features = dot.genes.malig, group.by = "ICA_mes",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  theme_bw() + scale_size(range = c(0, 6)) +
  theme(axis.text = element_text(size = 12)) + #, axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())+
  scale_colour_viridis_c(option = "magma", direction=-1) + labs(y="",x="")+ scale_size(range = c(0, 6), limits = c(25, 100)) + 
p2 = DotPlot(oligo, features = dot.genes.malig, group.by = "ICA_mes",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  theme_bw() + scale_size(range = c(0, 6)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_colour_viridis_c(option = "magma", direction=-1) + labs(y="",x="") + theme(axis.text=element_text(size=12)) + scale_size(range = c(0, 6), limits = c(25, 100)) + 
p1 + p2 + plot_layout(ncol = 2)
ggsave(paste0(fig.dir,"/dotplot_all_malig_cells_mes_by_cancer_not_scaled.pdf"), width = 18, height = 3.8)



#also make separate ones per cancer
p1 = DotPlot(myeloid.astro, features = dot.genes.mac.final, group.by = "cluster_annotation", col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  scale_colour_viridis_c(option = "magma", direction=-1) + labs(y="",x="") + theme_bw() + scale_size(range = c(0, 6), limits = c(25, 100)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text=element_text(size=12))
p2 = DotPlot(myeloid.oligo, features = dot.genes.mac.final, group.by = "cluster_annotation", col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  scale_colour_viridis_c(option = "magma", direction=-1) + labs(y="",x="") + theme_bw() + scale_size(range = c(0, 6), limits = c(25, 100)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text=element_text(size=12))
p1+p2
ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_by_cancer_scaled.pdf"), width = 18, height = 3.8)

# #combined dot plot with key LR
# DotPlot(myeloid.astro, features = dot.genes.mac.final, group.by = "grade",
#         col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
#   scale_colour_viridis_c(option = "magma", direction=-1) + labs(y="Grade",x="") + theme_bw() + scale_size(range = c(0, 6), limits = c(25, 100)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text=element_text(size=12))
# ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_combined_by_grade_scaled.pdf"), width = 8, height = 3.8)


oligo$ICA_score_state_grade = paste(oligo$ICA_score_states,"grade",oligo$grade)
astro$ICA_score_state_grade = paste(astro$ICA_score_states,"grade",astro$grade)
astro$ICA_score_state_grade = factor(astro$ICA_score_state_grade, levels = sort(unique(astro$ICA_score_state_grade)))
oligo$ICA_score_state_grade = factor(oligo$ICA_score_state_grade, levels = sort(unique(oligo$ICA_score_state_grade)))


Idents(astro) <- "ICA_score_state_grade"
Idents(oligo) <- "ICA_score_state_grade"
p1 = DotPlot(astro, features = dot.genes.malig, group.by = "ICA_score_state_grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p2 = DotPlot(oligo, features = dot.genes.malig, group.by = "ICA_score_state_grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p1 + p2 #+ plot_layout(guides = "collect")
ggsave(paste0(fig.dir,"/dotplot_all_malig_cells_by_cancer_by_ICA_score_state_by_grade_scaled.pdf"), width = 18, height = 3.5)

p1 = DotPlot(astro, features = dot.genes.malig, group.by = "ICA_score_state_grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p2 = DotPlot(oligo, features = dot.genes.malig, group.by = "ICA_score_state_grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p1 + p2 #+ plot_layout(guides = "collect")
ggsave(paste0(fig.dir,"/dotplot_all_malig_cells_by_cancer_by_ICA_score_state_by_grade_not_scaled.pdf"), width = 18, height = 3.5)




oligo$ICA_mes_grade = paste(oligo$ICA_mes,"grade",oligo$grade)
astro$ICA_mes_grade = paste(astro$ICA_mes,"grade",astro$grade)
astro$ICA_mes_grade = factor(astro$ICA_mes_grade, levels = sort(unique(astro$ICA_mes_grade)))
oligo$ICA_mes_grade = factor(oligo$ICA_mes_grade, levels = sort(unique(oligo$ICA_mes_grade)))


Idents(astro) <- "ICA_mes_grade"
Idents(oligo) <- "ICA_mes_grade"
p1 = DotPlot(astro, features = dot.genes.malig, group.by = "ICA_mes_grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p2 = DotPlot(oligo, features = dot.genes.malig, group.by = "ICA_mes_grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p1 + p2 #+ plot_layout(guides = "collect")
ggsave(paste0(fig.dir,"/dotplot_all_malig_cells_mes_by_cancer_by_grade_scaled.pdf"), width = 18, height = 4)

p1 = DotPlot(astro, features = dot.genes.malig, group.by = "ICA_mes_grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p2 = DotPlot(oligo, features = dot.genes.malig, group.by = "ICA_mes_grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p1 + p2 #+ plot_layout(guides = "collect")
ggsave(paste0(fig.dir,"/dotplot_all_malig_cells_mes_by_cancer_by_grade_not_scaled.pdf"), width = 18, height = 4)





#statistical test to show enrichment in grade 4 astrocytomas
df = all.myeloid.filt@meta.data[,c("orig.ident","disease","grade","cluster_annotation")]
colnames(df) = c("sample","subtype","grade","type")

#compute proportions of SPP1+ in each sample
df.prop <- df %>%
  group_by(sample, subtype, grade) %>%
  summarise(
    n_total = n(),
    n_target = sum(type == "Macrophage SPP1+"),
    prop = n_target / n_total,
    .groups = "drop"
  )

#assign to groups, either grade 4 astro or other
df.prop$group <- with(
  df.prop,
  ifelse(subtype == "astro" & grade == 4,
         "Astro grade 4",
         "Other")
)
df.prop$group = factor(df.prop$group, levels = c("Other","Astro grade 4"))


fit <- glmmTMB(
  cbind(n_target, n_total - n_target) ~ group,
  family = betabinomial(),
  data = df.prop
)
summary(fit)

#plot per grade and subtype
ggplot(df.prop, aes(x = interaction(subtype, grade), y = prop)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.02) +
  theme_classic() + ylab("Proportion of SPP1+ macrophages") + xlab("Subtype.grade")
ggsave(paste0(fig.dir,"/boxplot_proportion_SPP1_mac_per_grade_subtype.pdf"), width = 6, height = 4)





#figure to show what mes cells were originally defined as
orig.state.oligo = oligo$ICA_score_states[which(oligo$ICA_mes == "reactive")]
orig.state.astro = astro$ICA_score_states[which(astro$ICA_mes == "reactive")]


df = astro@meta.data[,c("ICA_score_states","ICA_mes")]
colnames(df) = c("original_state","new_state")

# All cells
all_bar <- df %>%
  count(original_state) %>%
  mutate(bar = "All cells", prop = n / sum(n))

# Cells that were redefined
redefined_bar <- df %>%
  filter(new_state == "reactive") %>%
  count(original_state) %>%
  mutate(bar = "Reactive/MES", prop = n / nrow(df))

plot_df <- bind_rows(all_bar, redefined_bar)

state_cols = gg_color_hue(4)
names(state_cols) = c("astro-like","intermediate","neuro-like","OPC-like")


p1 = ggplot(plot_df, aes(x = bar, y = prop, fill = original_state)) +
  geom_col(width = 0.65) + scale_fill_manual(values = state_cols) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  theme_classic() + guides(fill="none") + labs(x = NULL, y = "Proportion", fill = "Original state", title = "Astrocytoma")


df = oligo@meta.data[,c("ICA_score_states","ICA_mes")]
colnames(df) = c("original_state","new_state")

# All cells
all_bar <- df %>%
  count(original_state) %>%
  mutate(bar = "All cells", prop = n / sum(n))

# Cells that were redefined
redefined_bar <- df %>%
  filter(new_state == "reactive") %>%
  count(original_state) %>%
  mutate(bar = "Reactive/MES", prop = n / nrow(df))

plot_df <- bind_rows(all_bar, redefined_bar)

p2 = ggplot(plot_df, aes(x = bar, y = prop, fill = original_state)) +
  geom_col(width = 0.65) + scale_fill_manual(values = state_cols) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  theme_classic() + labs(x = NULL, y = "Proportion", fill = "Original state", title = "Oligodendroglioma")
p1+p2
ggsave(paste0(fig.dir,"/barplot_cell_states_redefined_as_mes.pdf"), width = 6, height = 4)










wnt.genes = c("LRP5","LRP6","FZD1","FZD3","FZD7",
              "CTNNB1","AXIN2","LEF1","TCF7",
              "ROR1","ROR2","RYK","FZD2","FZD5",
              "VANGL1","VANGL2","PRICKLE1")

p1 = DotPlot(astro, features = wnt.genes, group.by = "ICA_score_states",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p2 = DotPlot(oligo, features = wnt.genes, group.by = "ICA_score_states",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p1 + p2 #+ plot_layout(guides = "collect")
ggsave(paste0(fig.dir,"/dotplot_wnt_genes_all_malig_cells_by_cancer_by_ICA_score_state.pdf"), width = 20, height = 3.5)



Idents(astro) <- "ICA_score_state_grade"
Idents(oligo) <- "ICA_score_state_grade"


p1 = DotPlot(astro, features = wnt.genes, group.by = "ICA_score_state_grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p2 = DotPlot(oligo, features = wnt.genes, group.by = "ICA_score_state_grade",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=F) +
  RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
  scale_colour_viridis_c(option = "magma", direction=-1)
p1 + p2 #+ plot_layout(guides = "collect")
ggsave(paste0(fig.dir,"/dotplot_wnt_genes_all_malig_cells_by_cancer_by_ICA_score_state_by_grade.pdf"), width = 22, height = 3.5)






####identify regulons for macrophage genes
#first in astro
regulons = as.data.frame(fread("/mnt/beegfs/home/alyne/stageout/scenicplus_op_newdb_filt_astro/outs/eRegulon_direct.tsv",
                 data.table=F))

out.dir = paste0("/mnt/beegfs/home/alyne/stageout/scenicplus_op_newdb_filt_astro/figures_scenicplus_newdb_astro")
use.extended = F
#read in region AUC scores
TF.AUCg.d = fread(paste0(out.dir,"/direct_gene_based_AUC_astro.tsv.gz"), header=T, data.table = F)
if(use.extended){
  TF.AUCg.e = fread(paste0(out.dir,"/extended_gene_based_AUC_astro.tsv.gz"), header=T, data.table = F)
}else{
  TF.AUCg.e = data.frame(Cell = TF.AUCg.d$Cell)
}
TF.AUCg = cbind(TF.AUCg.d,TF.AUCg.e[,-1])
rownames(TF.AUCg) = TF.AUCg$Cell
rownames(TF.AUCg) = sapply(rownames(TF.AUCg), function(x) {unlist(strsplit(x,"___",fixed=T))[1]})
TF.AUCg$Cell = NULL
# rownames(TF.AUCg) = sapply(rownames(TF.AUCg), function(x) {
#   parts = unlist(strsplit(x,"-",fixed=T))
#   paste(parts[3],paste0(parts[1],"-",parts[2]),sep="#")
# })
# TF.AUCg.ss = TF.AUCg[names.both,]

## gene expression
expr.mat <- GetAssayData(all.myeloid.filt, assay = "SCT", layer = "data")

## make sure cells are in the same order
common.cells <- intersect(colnames(expr.mat), rownames(TF.AUCg))

expr.mat <- expr.mat[, common.cells]
TF.AUCg.ss <- TF.AUCg[common.cells, ]

#only keep +/+
TF.AUCg.ss = TF.AUCg.ss[,grepl("+/+",colnames(TF.AUCg.ss),fixed=T)]


#grade 4 mac population
genes.of.int = c("SPP1","PLTP","ACSL1","MAP3K8","VSIG4")

#for each gene get whcih TF regulons it is in
TF.list = sapply(genes.of.int, function(x) {
  unique(regulons$TF[which(regulons$Gene==x)])
})

all.tfs <- sort(unique(unlist(TF.list)))
all.tfs <- all.tfs[which(!(all.tfs %in% c("MAF","ATF6","FOSL2","ID4","IRF8","IKZF1","KLF6","MAFB","RUNX2")))] #remove MAF because it's not in scenicplus output

plot.df <- expand.grid(Gene = genes.of.int, TF = all.tfs, stringsAsFactors = FALSE)
plot.df$mean_rank <- NA_real_
plot.df$n_regions <- NA_integer_

for(i in seq_len(nrow(plot.df))){
  
  g <- plot.df$Gene[i]
  tf <- plot.df$TF[i]
  
  rel <- regulons[regulons$Gene == g & regulons$TF == tf, ]
  
  if(nrow(rel) > 0){
    plot.df$mean_rank[i] <- mean(sort(rel$triplet_rank)[1:min(2,nrow(rel))])
    plot.df$n_regions[i] <- nrow(rel)
  }
}

plot.df$cor_activity   <- NA_real_
plot.df$cor_expression <- NA_real_

for(i in seq_len(nrow(plot.df))){
  g  <- plot.df$Gene[i]
  tf <- plot.df$TF[i]
  tf.col.name = colnames(TF.AUCg.ss)[grepl(paste0(tf,"_direct"),colnames(TF.AUCg.ss))]
  if(tf %in% TF.list[[g]] & g %in% rownames(expr.mat) & tf.col.name %in% colnames(TF.AUCg.ss)){
    gene.expr <- as.numeric(expr.mat[g, ])
    ## TF activity
    tf.act <- TF.AUCg.ss[, tf.col.name]
    plot.df$cor_activity[i] <-
      suppressWarnings(cor(tf.act, gene.expr, method = "spearman", use = "pairwise.complete.obs"))
    ## TF expression
    if(tf %in% rownames(expr.mat)){
      tf.expr <- as.numeric(expr.mat[tf, ])
      
      plot.df$cor_expression[i] <-
        suppressWarnings(cor(tf.expr, gene.expr, method = "spearman", use = "pairwise.complete.obs"))
    }
  }
}

plot.df$Gene <- factor(plot.df$Gene, levels = rev(genes.of.int))
plot.df$TF <- factor(plot.df$TF, levels = all.tfs)

plot.df$cor_expression[which(plot.df$cor_expression>0.99999)] = NA_real_

p1 <- ggplot(plot.df, aes(TF, Gene, fill = mean_rank)) +
  geom_tile(colour = "grey90") + scale_fill_viridis_c(option = "magma", na.value = "white", trans = "log10",) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_blank()) +
  labs(fill = "Mean\ntriplet rank") + theme(text=element_text(size=20)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p2 <- ggplot(plot.df, aes(TF, Gene, fill = n_regions)) +
  geom_tile(colour = "grey90") + scale_fill_viridis_c(option = "plasma", na.value = "white", direction = -1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_blank()) +
  labs(fill = "# regions") + theme(text=element_text(size=20)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p3 <- ggplot(plot.df, aes(TF, Gene, fill = cor_activity)) +
  geom_tile(colour = "grey90") + scale_fill_viridis_c(option = "plasma", na.value = "white", direction = -1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90,     hjust = 0.5,
                                                vjust = 0.5), panel.grid = element_blank()) +
  labs(fill = "Act. correlation") + theme(text=element_text(size=20))
# p4 <- ggplot(plot.df, aes(TF, Gene, fill = cor_expression)) +
#   geom_tile(colour = "grey90") + scale_fill_viridis_c(option = "plasma", na.value = "white", direction = -1) +
#   theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_blank()) +
#   labs(fill = "Expr. correlation") + theme(text=element_text(size=20))
wrap_plots(list(p1,p2,p3), ncol=1)
ggsave(paste0(fig.dir,"/heatmap_spp1_TFs_astro.pdf"), width=8, height=7)



#



#TFs of interest
TFs.of.int = unique(c("CEBPB","CEBPD","FLI1",
               "MITF","FLI1","ETS2","CEBPB",
               "FLI1","RBPJ","MAFB","ETV6",
               "ELF2","ETS2","RUNX1","CEBPB"))





#then in oligo
regulons = as.data.frame(fread("/Users/alyne/Library/CloudStorage/Dropbox/Curie/Curie_StCloud/anne-marie/10x_RNA_ATAC_oligo_analysis/SCENICplus_soroban/oligo_full_run/outs_op_newdb_filt/eRegulon_direct_oligo.tsv",
                               data.table=F))

genes.of.int = c("CD163","SPP1","MRC1","PLTP","ACSL1","C1QB","NAMPT","SLC11A1")

#for each gene get whcih TF regulons it is in
TF.list = sapply(genes.of.int, function(x) {
  unique(regulons$TF[which(regulons$Gene==x)])
})

#for each TF, if has three or more regions, find average rank of top 3
lapply(seq(length(TF.list)), function(x) {
  if(length(TF.list[[x]])>0){
    sort(sapply(TF.list[[x]], function(y) {
      rel.rows = subset(regulons, Gene == genes.of.int[x] & TF == y)
      if(nrow(rel.rows)>=1){
        mean(sort(rel.rows$triplet_rank)[1:1])
      }else{
        NA_real_
      }
    }),na.last=T)
  }else{
    c(NA)
  }
})






library(ArchR)
library(parallel)


cancer = "astro"
base.dir = "/mnt/beegfs/home/alyne"
setwd(base.dir)
out.dir = paste0("stageout/ArchR_output_",cancer)
dir.create(out.dir)
N.min.grp = 100

set.seed(1)

addArchRLocking(locking = T)
addArchRThreads(threads = 1)
addArchRGenome("hg38")

#read in sample list
samps = fread(paste0("stagein/scenicplus_",cancer,"/sample_list.txt"),data.table = F, header = F)[,1]
samp.md = as.data.frame(read_excel("~/persistent/sample_summary_astro_oligo_anon.xlsx"))

oligo.ArchR = loadArchRProject(path = paste0(out.dir,"/project1"), force = T, showLogo = F)

#read in dataframe with trajectory 
df.traj = fread(paste0("persistent/trajectory/df_file_traj_",cancer,".tsv.gz"), data.table = F)
#add archr style names to df.traj
df.traj$cell.names.archr = sapply(seq(nrow(df.traj)), function(x) {
  cell.n = paste(unlist(strsplit(df.traj$cell[x],"-"))[1:2],collapse="-")
  paste(df.traj$sample[x],cell.n,sep="#")
})

anno = fread(paste0(fig.dir,"/cell_annotations_minus_S0009.tsv.gz"), data.table=F)
#get names that match ArchR
anno$cell.name.archr <- sub("^(.+)-1-(.+)$","\\2#\\1-1",anno$cell.name)

myeloid.cells <- intersect(anno$cell.name.archr,getCellNames(oligo.ArchR))

# oligo.myeloid <- subsetArchRProject(
#   ArchRProj = oligo.ArchR,
#   cells = myeloid.cells,
#   outputDirectory = "oligo_myeloid",
#   force = TRUE
# )

oligo.myeloid <- oligo.ArchR[myeloid.cells, ]

oligo.myeloid <- addBgdPeaks(ArchRProj = oligo.myeloid, force=T)
oligo.myeloid <- addDeviationsMatrix(ArchRProj = oligo.myeloid, peakAnnotation = "Motif", force = TRUE)
motifs <- getMatrixFromProject(ArchRProj = oligo.myeloid, useMatrix = "MotifMatrix")
rownames(motifs)[grep("FLI1|CEBPB|CEBPD|IKZF1|MAFB|ELF1|RUNX1|PPARG",rownames(motifs))]
TFs <- c("FLI1","CEBPB","CEBPD","IKZF1","MAFB","ELF1","RUNX1","PPARG")

motif.mat <- motifs[TFs,]
anno2 <- anno[match(colnames(motif.mat),anno$cell.name.archr),]

#compute mean accessiblity per TF
library(Matrix)

motif.df <- as.data.frame(t(as.matrix(motif.mat)))

motif.df$cluster <- anno2$cluster

motif.mean <- aggregate(
  . ~ cluster,
  data = motif.df,
  FUN = mean
)

rownames(motif.mean) <- motif.mean$cluster
motif.mean$cluster <- NULL

motif.mean <- as.matrix(motif.mean)

#heatmap
Heatmap(
  t(scale(t(motif.mean))),
  name = "chromVAR\n(z-score)",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left"
)

###########################
## process to get markers for strange samples
###########################

#astro_P0027_ini_S0052
#just remove the samples with IFN or hypoxic
samps.to.rem = c("astro_P0026_ini_S0050", "oligo_P0004_ini_S0009")
all.myeloid.filt.v2 = subset(all.myeloid, subset = !(orig.ident %in% samps.to.rem))

all.myeloid.filt.v2 <- IntegrateLayers(
  object = all.myeloid.filt.v2, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.har",
  verbose = FALSE, assay = "RNA"
)
all.myeloid.filt.v2 <- FindNeighbors(all.myeloid.filt.v2, reduction = "integrated.har", dims = 1:30, graph.name = "har_nn")
all.myeloid.filt.v2 <- RunUMAP(all.myeloid.filt.v2, reduction = "integrated.har", dims = 1:30, reduction.name = "umap.har")
res = 0.3
all.myeloid.filt.v2 <- FindClusters(all.myeloid.filt.v2, resolution = res, cluster.name = paste0("har_clusters",res), graph = "har_nn")
DimPlot(all.myeloid.filt.v2, group.by=c(paste0("har_clusters",res),"disease","myeloid_CellType","label_program","grade","orig.ident"), 
        reduction = "umap.har", shuffle = T)
ggsave(paste0(fig.dir,"/umaps_filt_myeloid_cells_clusters_with_S0052.pdf"), width=25, height=8)

FeaturePlot(all.myeloid.filt.v2, features = c("Microglia","Complement_Immunosuppressive","Neutrophil","Scavenger_Immunosuppressive","G2_M","IFN_Response",
                                      "Inflammatory_microglia","Macrophage","IL1B_CD83_TNF_Inflammatory","cDC","G1_S","Hypoxia","HS_UPR","Monocyte"),
            reduction = "umap.har")
ggsave(paste0(fig.dir,"/umaps_filt_myeloid_cells_miller_scores_with_S0052.pdf"), width=25, height=25)


#get markers of cluster 3 (S0052) compared to cluster 2 (which is also immunoregulatory/suppressive)
all.myeloid.filt.v2 = PrepSCTFindMarkers(all.myeloid.filt.v2)
markers.52 = FindMarkers(all.myeloid.filt.v2, group.by = paste0("har_clusters",res),ident.1 = 3, ident.2 = 2)
genes52.to.plot = rownames(head(markers.52[markers.52$avg_log2FC>0,],N.per.cat))

#remove genes already in main cluster marker genes
genes52.to.plot = genes52.to.plot[!(genes52.to.plot %in% genes.to.plot)]
cells.52 = names(all.myeloid.filt.v2$har_clusters0.3[which(all.myeloid.filt.v2$har_clusters0.3==3)])
cells.52.to.plot = cells.52[sample(length(cells.52),min(Ncell,length(cells.52)))]
cells.52.to.plot = cells.52.to.plot











#astro_P0026_ini_S0050
#sample S0050 has two separate clusters, so do clustering on Miller scores to identify cells in each cluster
#then use full object to get marker genes for each cluster
my.scores = all.myeloid@meta.data[all.myeloid$disease=="astro",c("Microglia","Complement_Immunosuppressive","Neutrophil","Scavenger_Immunosuppressive","G2_M","IFN_Response",
                                                                 "Inflammatory_microglia","Macrophage","IL1B_CD83_TNF_Inflammatory","cDC","G1_S","Hypoxia","HS_UPR","Monocyte")]
assay.v5 <- CreateAssay5Object(data = scale(t(my.scores)))
all.myeloid_slim <- CreateSeuratObject(assay.v5)
all.genes = rownames(all.myeloid_slim)
#also dim reduction
colnames(my.scores) = paste0("SIXD_",1:ncol(my.scores))
all.myeloid_slim[["sixD"]] <- CreateDimReducObject(embeddings = as.matrix(my.scores),key = "SIXD_", assay = "RNA")
all.myeloid_slim <- FindNeighbors(all.myeloid_slim, reduction = "sixD",
                                  dims = 1:ncol(my.scores), graph.name = c("ICA_state_graph_nn","ICA_state_graph_snn"))
#cluster
res = 0.45
all.myeloid_slim <- FindClusters(all.myeloid_slim, resolution = res, cluster.name = paste0("ICA_state_clusters",res), 
                                 graph = "ICA_state_graph_nn", group.singletons = TRUE)
all.myeloid_slim$orig.ident = all.myeloid$orig.ident[which(all.myeloid$disease=="astro")]

all.myeloid.astro = subset(all.myeloid, subset = disease=="astro")
all.myeloid.astro$Miller_clusters = all.myeloid_slim$ICA_state_clusters0.45

DefaultAssay(all.myeloid.astro) <- "RNA"
all.myeloid.astro <- NormalizeData(all.myeloid.astro)
all.myeloid.astro <- FindVariableFeatures(all.myeloid.astro)
all.myeloid.astro <- JoinLayers(all.myeloid.astro)
markers.astro <- FindAllMarkers(all.myeloid.astro, assay = "RNA", group.by = "Miller_clusters")

#get top genes for the two relevant clusters
genes.to.plot.50 = unlist(lapply(c(4,5), function(x) {
  markers.astro[markers.astro$cluster==x & markers.astro$avg_log2FC>0,"gene"][1:N.per.cat]
}))
genes.to.plot.50 = genes.to.plot.50[!(genes.to.plot.50 %in% c(genes52.to.plot,genes.to.plot,genes09.to.plot))]
cells.to.plot.50.list = lapply(c(4,5), function(x) {
  all.cells.cl = names(all.myeloid.astro$Miller_clusters[which(all.myeloid.astro$Miller_clusters==x & 
                                                                 all.myeloid.astro$orig.ident == "astro_P0026_ini_S0050")])
  all.cells.cl[sample(length(all.cells.cl),min(Ncell,length(all.cells.cl)))]
})
cells.to.plot.50 = unlist(cells.to.plot.50.list)




mat = all.myeloid@assays$SCT@data
mat.scaled <- t(scale(t(mat)))
all.cells.plot = c(cells.to.plot,cells.52.to.plot,cells.to.plot.50,cells.09.to.plot)
all.genes.plot = c(genes.to.plot,genes52.to.plot,genes.to.plot.50,genes09.to.plot)
mat.to.plot = mat.scaled[all.genes.plot,all.cells.plot]

col.annot.df = data.frame(cluster = c(rep("Microglia",Ncell), rep("Inflammatory",Ncell), rep("Immunoregulatory",Ncell), rep("Cycling",length(cells.to.plot.list[[4]])),
                                      rep("Complement",length(cells.52.to.plot)), rep("Immunoregulatory (oxidative)",length(cells.to.plot.50.list[[1]])), 
                                      rep("Hypoxic",length(cells.to.plot.50.list[[2]])),  rep("IFN-responsive",length(cells.09.to.plot))))
col.annot.df$cluster = factor(col.annot.df$cluster, levels = c("Microglia","Inflammatory","Immunoregulatory", "Cycling",
                                                               "Complement","Immunoregulatory (oxidative)",
                                                               "Hypoxic","IFN-responsive"))
rownames(col.annot.df) = all.cells.plot
cols.cluster = gg_color_hue(7)
names(cols.cluster) = levels(col.annot.df$cluster)
cols.list.col = list(cluster = cols.cluster)
top_anno <- HeatmapAnnotation(df = col.annot.df, col = cols.list.col,
                              annotation_legend_param = list(cluster = list(labels_gp = gpar(fontsize = 20),
                                                                            title_gp  = gpar(fontsize = 22))))

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht <- Heatmap(
  mat.to.plot,
  col = col_fun,
  name = "Gene expression",
  cluster_rows = F,
  cluster_columns = F,
  top_annotation  = top_anno,
  show_row_names = T,
  show_column_names = F,
  row_title = NULL,
  column_title = NULL,
  show_column_dend = F,
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect( x, y, w + unit(0.3, "pt"), h + unit(0.3, "pt"),
               gp = gpar(fill = fill, col = NA) )
  },  heatmap_legend_param = list(legend_height = unit(6, "cm"), labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 22))
)
pdf(paste0(fig.dir,"/heatmap_gex_myeloid_full_nocluster.pdf"), width = 20, height = 20, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()



#do clustered version just to see how it looks
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht <- Heatmap(
  mat.to.plot,
  col = col_fun,
  name = "Gene expression",
  cluster_rows = T,
  cluster_columns = T,
  top_annotation  = top_anno,
  show_row_names = T,
  show_column_names = F,
  row_title = NULL,
  column_title = NULL,
  show_column_dend = F,
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect( x, y, w + unit(0.3, "pt"), h + unit(0.3, "pt"),
               gp = gpar(fill = fill, col = NA) )
  },
  heatmap_legend_param = list(legend_height = unit(6, "cm"), labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 22))
)
pdf(paste0(fig.dir,"/heatmap_gex_myeloid_full_clustered.pdf"), width = 20, height = 20, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


























# 
# #S0050 and S0009 defeinitely different and need to be removed
# #S0052 overlaps with main cluster 2, so need to find DE genes relative to this and use to see how different
# 
# 
# 
# 
# 
# 
# 
# 
# lapply(c(0:3), function(x) {
#   head(markers[markers$cluster==x & markers$avg_log2FC>0,],30)
# })
# #3 = cell cycle, 2 = Resident-like tissue-remodelling TAMs, 1 = activated inflammatory, 0 = homeostatic microglia
# 
# 
# 
# #plot heatmap with top 20 genes from each state
# #also include top 20 genes from each of the states in specific samples
# genes.to.plot = unlist(lapply(c(0:3), function(x) {
#   markers[markers$cluster==x & markers$avg_log2FC>0,"gene"][1:30]
# }))
# 
# #subsample 1000 cells per cluster to plot, plus all cycling cells
# cells.to.plot = unlist(lapply(c(0:2), function(x) {
#   all.cells.cl = names(all.myeloid.filt$har_clusters0.3[which(all.myeloid.filt$har_clusters0.3==x)])
#   all.cells.cl[sample(length(all.cells.cl),1000)]
# }))
# cells.to.plot = c(cells.to.plot, 
#                   names(all.myeloid.filt$har_clusters0.3[which(all.myeloid.filt$har_clusters0.3==3)]))
# mat = all.myeloid.filt@assays$SCT@data
# mat.scaled <- t(scale(t(mat)))
# mat.to.plot = mat.scaled[genes.to.plot,cells.to.plot]
# 
# col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
# ht <- Heatmap(
#   mat.to.plot,
#   col = col_fun,
#   name = "Gene expression",
#   cluster_rows = F,
#   cluster_columns = F,
#   # row_split = split_groups.filt,       # vector of cluster IDs
#   # row_gap = unit(0.5, "mm"),
#   # column_gap = unit(0.5, "mm"),
#   show_row_names = T,
#   show_column_names = F,
#   row_title = NULL,
#   column_title = NULL,
#   show_column_dend = F,
#   layer_fun = function(j, i, x, y, w, h, fill) {
#     grid.rect( x, y, w + unit(0.3, "pt"), h + unit(0.3, "pt"),
#                gp = gpar(fill = fill, col = NA) )
#   }
# )
# pdf(paste0(fig.dir,"/heatmap_gex_myeloid_filt_nocluster.pdf"), width = 15, height = 15, useDingbats = FALSE)
# draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
# dev.off()
# 
# 
# 
# #now need to check specifically S0009 in oligo, do all cells end in one cluster, 
# #and if so, which genes characterize that cluster
# #make new assay with relevant scores
# my.scores = all.myeloid@meta.data[all.myeloid$disease=="oligo",c("Microglia","Complement_Immunosuppressive","Neutrophil","Scavenger_Immunosuppressive","G2_M","IFN_Response",
#                                                                  "Inflammatory_microglia","Macrophage","IL1B_CD83_TNF_Inflammatory","cDC","G1_S","Hypoxia","HS_UPR","Monocyte")]
# assay.v5 <- CreateAssay5Object(data = scale(t(my.scores)))
# all.myeloid_slim <- CreateSeuratObject(assay.v5)
# all.genes = rownames(all.myeloid_slim)
# #also dim reduction
# colnames(my.scores) = paste0("SIXD_",1:ncol(my.scores))
# all.myeloid_slim[["sixD"]] <- CreateDimReducObject(
#   embeddings = as.matrix(my.scores),
#   key = "SIXD_",
#   assay = "RNA"
# )
# all.myeloid_slim <- FindNeighbors(all.myeloid_slim, reduction = "sixD",
#                                   dims = 1:ncol(my.scores), graph.name = c("ICA_state_graph_nn","ICA_state_graph_snn")
# )
# #cluster
# res = 0.4
# all.myeloid_slim <- FindClusters(all.myeloid_slim, resolution = res, cluster.name = paste0("ICA_state_clusters",res), 
#                                  graph = "ICA_state_graph_nn", group.singletons = TRUE)
# all.myeloid_slim$orig.ident = all.myeloid$orig.ident[which(all.myeloid$disease=="oligo")]
# table(all.myeloid_slim$ICA_state_clusters0.4[which(all.myeloid_slim$orig.ident=="oligo_P0004_ini_S0009")])
# #2014/2545 cells in cluster 2
# 
# 
# all.myeloid.oligo = subset(all.myeloid, subset = disease=="oligo")
# all.myeloid.oligo$Miller_clusters = all.myeloid_slim$ICA_state_clusters0.4
# all.myeloid.oligo = PrepSCTFindMarkers(all.myeloid.oligo)
# markers.oligo = FindAllMarkers(all.myeloid.oligo, group.by = "Miller_clusters")
# 
# S0009.markers = markers.oligo[markers.oligo$cluster==2 & markers.oligo$avg_log2FC>0,"gene"][1:30]
# cells.to.plot.S0009 = unlist(lapply(c(2), function(x) {
#   all.cells.cl = names(all.myeloid.oligo$Miller_clusters[which(all.myeloid.oligo$Miller_clusters==x & 
#                                                                  all.myeloid.oligo$orig.ident == "oligo_P0004_ini_S0009")])
#   all.cells.cl[sample(length(all.cells.cl),min(1000,length(all.cells.cl)))]
# }))
# 
# #now do the same for astro samples
# my.scores = all.myeloid@meta.data[all.myeloid$disease=="astro",c("Microglia","Complement_Immunosuppressive","Neutrophil","Scavenger_Immunosuppressive","G2_M","IFN_Response",
#                                                                  "Inflammatory_microglia","Macrophage","IL1B_CD83_TNF_Inflammatory","cDC","G1_S","Hypoxia","HS_UPR","Monocyte")]
# assay.v5 <- CreateAssay5Object(data = scale(t(my.scores)))
# all.myeloid_slim <- CreateSeuratObject(assay.v5)
# all.genes = rownames(all.myeloid_slim)
# #also dim reduction
# colnames(my.scores) = paste0("SIXD_",1:ncol(my.scores))
# all.myeloid_slim[["sixD"]] <- CreateDimReducObject(
#   embeddings = as.matrix(my.scores),
#   key = "SIXD_",
#   assay = "RNA"
# )
# all.myeloid_slim <- FindNeighbors(all.myeloid_slim, reduction = "sixD",
#                                   dims = 1:ncol(my.scores), graph.name = c("ICA_state_graph_nn","ICA_state_graph_snn")
# )
# #cluster
# res = 0.45
# all.myeloid_slim <- FindClusters(all.myeloid_slim, resolution = res, cluster.name = paste0("ICA_state_clusters",res), 
#                                  graph = "ICA_state_graph_nn", group.singletons = TRUE)
# all.myeloid_slim$orig.ident = all.myeloid$orig.ident[which(all.myeloid$disease=="astro")]
# table(all.myeloid_slim$ICA_state_clusters0.45[which(all.myeloid_slim$orig.ident=="astro_P0026_ini_S0050")])
# #543 in cluster 4 and 211 in cluster 5 out of 866
# table(all.myeloid_slim$ICA_state_clusters0.45[which(all.myeloid_slim$orig.ident=="astro_P0027_ini_S0052")])
# #1302 in cluster 1 out of 1360
# 
# all.myeloid.astro = subset(all.myeloid, subset = disease=="astro")
# all.myeloid.astro$Miller_clusters = all.myeloid_slim$ICA_state_clusters0.45
# # all.myeloid.astro = PrepSCTFindMarkers(all.myeloid.astro)
# # markers.astro = FindAllMarkers(all.myeloid.astro, group.by = "Miller_clusters")
# 
# 
# DefaultAssay(all.myeloid.astro) <- "RNA"
# all.myeloid.astro <- NormalizeData(all.myeloid.astro)
# all.myeloid.astro <- FindVariableFeatures(all.myeloid.astro)
# all.myeloid.astro <- JoinLayers(all.myeloid.astro)
# markers.astro <- FindAllMarkers(all.myeloid.astro, assay = "RNA", group.by = "Miller_clusters")
# 
# genes.to.plot.astro = unlist(lapply(c(1,4,5), function(x) {
#   markers.astro[markers.astro$cluster==x & markers.astro$avg_log2FC>0,"gene"][1:30]
# }))
# 
# cells.to.plot.astro = unlist(lapply(c(1,4,5), function(x) {
#   all.cells.cl = names(all.myeloid.astro$Miller_clusters[which(all.myeloid.astro$Miller_clusters==x & 
#                                                                  all.myeloid.astro$orig.ident %in% c("astro_P0026_ini_S0050","astro_P0027_ini_S0052"))])
#   all.cells.cl[sample(length(all.cells.cl),min(1000,length(all.cells.cl)))]
# }))
# 
# 
# 
# mat = all.myeloid@assays$SCT@data
# mat.scaled <- t(scale(t(mat)))
# mat.to.plot = mat.scaled[c(genes.to.plot,S0009.markers,genes.to.plot.astro),
#                          c(cells.to.plot,cells.to.plot.S0009,cells.to.plot.astro)]
# 
# col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
# ht <- Heatmap(
#   mat.to.plot,
#   col = col_fun,
#   name = "Gene expression",
#   cluster_rows = F,
#   cluster_columns = F,
#   # row_split = split_groups.filt,       # vector of cluster IDs
#   # row_gap = unit(0.5, "mm"),
#   # column_gap = unit(0.5, "mm"),
#   show_row_names = T,
#   show_column_names = F,
#   row_title = NULL,
#   column_title = NULL,
#   show_column_dend = F,
#   layer_fun = function(j, i, x, y, w, h, fill) {
#     grid.rect( x, y, w + unit(0.3, "pt"), h + unit(0.3, "pt"),
#                gp = gpar(fill = fill, col = NA) )
#   }
# )
# pdf(paste0(fig.dir,"/heatmap_gex_myeloid_filt_and_specific_nocluster.pdf"), width = 20, height = 25, useDingbats = FALSE)
# draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
# dev.off()
# 
# 
# 
# #S0050 and S0009 defeinitely different and need to be removed
# #S0052 overlaps with main cluster 2, so need to find DE genes relative to this and use to see how different
# 
# 
# 
# 
# "astro_P0026_ini_S0050","astro_P0027_ini_S0052"
# 
# all.myeloid.oligo = subset(all.myeloid, subset = disease=="oligo")
# all.myeloid.oligo$Miller_clusters = all.myeloid_slim$ICA_state_clusters0.4
# all.myeloid.oligo = PrepSCTFindMarkers(all.myeloid.oligo)
# markers.oligo = FindAllMarkers(all.myeloid.oligo, group.by = "Miller_clusters")
# 
# S0009.markers = markers.oligo[markers.oligo$cluster==2 & markers.oligo$avg_log2FC>0,"gene"][1:30]
# 
# 
# 
# 
# 
# 
# 
# #clustering and UMAP
# all.myeloid <- FindNeighbors(all.myeloid, reduction = "integrated.har", dims = 1:50, graph.name = "har_nn")
# all.myeloid <- RunUMAP(all.myeloid, reduction = "integrated.har", dims = 1:50, reduction.name = "umap.har2")
# 
# res = 0.3
# all.myeloid <- FindClusters(all.myeloid, resolution = res, cluster.name = paste0("har_clusters",res), graph = "har_nn")
# 
# 
# DimPlot(all.myeloid, group.by=c(paste0("har_clusters",res),"disease","myeloid_CellType","label_program","grade","orig.ident"), 
#         reduction = "umap.har2", shuffle = T)
# ggsave(paste0(fig.dir,"/umaps_original_myeloid_cells_clusters.pdf"),
#        width=25, height=8)
# 
# FeaturePlot(all.myeloid, features = c("Microglia","Complement_Immunosuppressive","Neutrophil","Scavenger_Immunosuppressive","G2_M","IFN_Response",
#                                       "Inflammatory_microglia","Macrophage","IL1B_CD83_TNF_Inflammatory","cDC","G1_S","Hypoxia","HS_UPR","Monocyte"),
#             reduction = "umap.har2")
# ggsave(paste0(fig.dir,"/umaps_original_myeloid_cells_miller_scores.pdf"),
#        width=25, height=25)
# 
# #oligo_P0004_ini_S0009 and astro_P0026_ini_S0050  astro_P0027_ini_S0052
# 
# #make new assay with relevant scores
# my.scores = all.myeloid@meta.data[all.myeloid$disease=="astro",c("Microglia","Complement_Immunosuppressive","Neutrophil","Scavenger_Immunosuppressive","G2_M","IFN_Response",
#                          "Inflammatory_microglia","Macrophage","IL1B_CD83_TNF_Inflammatory","cDC","G1_S","Hypoxia","HS_UPR","Monocyte")]
# 
# 
# 
# 
# 
# 
# assay.v5 <- CreateAssay5Object(data = scale(t(my.scores)))
# all.myeloid_slim <- CreateSeuratObject(assay.v5)
# 
# all.genes = rownames(all.myeloid_slim)
# 
# #also dim reduction
# colnames(my.scores) = paste0("SIXD_",1:ncol(my.scores))
# all.myeloid_slim[["sixD"]] <- CreateDimReducObject(
#   embeddings = as.matrix(my.scores),
#   key = "SIXD_",
#   assay = "RNA"
# )
# 
# 
# all.myeloid_slim <- FindNeighbors(all.myeloid_slim, reduction = "sixD",
#   dims = 1:ncol(my.scores), graph.name = c("ICA_state_graph_nn","ICA_state_graph_snn")
# )
# 
# #cluster
# res = 0.3
# all.myeloid_slim <- FindClusters(all.myeloid_slim, resolution = res, cluster.name = paste0("ICA_state_clusters",res), 
#                          graph = "ICA_state_graph_nn", group.singletons = TRUE)
# #umap
# all.myeloid_slim = RunUMAP(all.myeloid_slim, features = all.genes,
#                    reduction.name = "UMAP_ICA_state", assay="RNA", slot="data")
# 
# FeaturePlot(all.myeloid_slim, features=all.genes, reduction = "UMAP_ICA_state")
# ggsave(paste0(fig.dir,"/umaps_myeloid_cells_all_Miller_astro.pdf"),width=25, height=25)
# 
# 
# #add annotations to plot
# all.myeloid_slim$disease = all.myeloid$disease[which(all.myeloid$disease=="astro")]
# all.myeloid_slim$myeloid_CellType = all.myeloid$myeloid_CellType
# all.myeloid_slim$label_program = all.myeloid$label_program
# all.myeloid_slim$grade = all.myeloid$grade
# all.myeloid_slim$orig.ident = all.myeloid$orig.ident
# 
# 
# DimPlot(all.myeloid_slim, group.by=c(paste0("ICA_state_clusters",res),"disease","myeloid_CellType","label_program","grade","orig.ident"), 
#         reduction = "UMAP_ICA_state", shuffle = T, label=T)
# ggsave(paste0(fig.dir,"/umaps_myeloid_cells_clusters.pdf"),
#        width=25, height=8)
# 
# 
# 
# #add Miller-based clusters back to full object
# all.myeloid$Miller_clusters = all.myeloid_slim$ICA_state_clusters0.3
# all.myeloid = PrepSCTFindMarkers(all.myeloid)
# markers = FindAllMarkers(all.myeloid, group.by = "Miller_clusters")
# 
# 
# all.myeloid[["UMAP_ICA_state"]] <- all.myeloid_slim[["UMAP_ICA_state"]]
# #microglia markers
# mic.mark = c("P2RY12","P2RY13","TMEM119","GPR34","SALL1","NAV3","FCHO2",
#              "AHR","SASH1")
# 
# 
# #macrophage markers
# mac.mark = c("FCGR3A","LILRB1","VCAN","S100A8","S100A9","CTSB",
#             "TYMP","SAT1")
# mac.genesets = list(mic.mark = mic.mark, mac.mark = mac.mark)
# 
# all.myeloid = AddModuleScore(all.myeloid, features = mac.genesets, name = names(mac.genesets))
# 
# FeaturePlot(all.myeloid, features=mic.mark, reduction = "UMAP_ICA_state")
# ggsave(paste0(fig.dir,"/umaps_myeloid_cells_all_mic_mark.pdf"),width=25, height=25)
# 
# FeaturePlot(all.myeloid, features=mac.mark, reduction = "UMAP_ICA_state")
# ggsave(paste0(fig.dir,"/umaps_myeloid_cells_all_mac_mark.pdf"),width=25, height=25)
# 
# FeaturePlot(all.myeloid, features=c("mic.mark1","mac.mark2"), reduction = "UMAP_ICA_state")
# ggsave(paste0(fig.dir,"/umaps_myeloid_cells_all_mark_genesets.pdf"),width=12, height=5)
# 
# 
# markers.filt = markers[markers$p_val_adj<0.01 & markers$avg_log2FC>0.8,]
# 
# 
# #list of ligands to check
# mac.ligs = c("TLR4","VSIR","SPP1",
#              "ITGA9","ITGB1","PTPRM","ITGB5","ITGAV","WNT5A",
#              "PDGFC","PDGFB","HLA-E","SLC1A3",
#              "SORL1","CD74","IL18","TLR2")
# 
# lapply(mac.ligs, function(x) {markers.filt[markers.filt$gene==x,]})
# 
# 
# #genes for dot plot
# dot.genes.mac = c("SPP1","ITGB1","ITGAV","ITGA9","ITGB5","PTPRM","WNT5A","TLR4","PDGFB","PDGFC","HLA-E","ABCA1")
# 
# 
# library(scCustomize)
# 
# myeloid.astro = subset(all.myeloid, subset = disease == "astro")
# myeloid.oligo = subset(all.myeloid, subset = disease == "oligo")
# 
# DotPlot(all.myeloid, features = dot.genes.mac, group.by = "Miller_clusters",
#   cols = c("grey95", "#b2182b")) +
#   RotatedAxis() + theme_bw()
# ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells.pdf"), width = 6, height = 2.7)
# 
# 
# Idents(myeloid.astro) <- "Miller_clusters"
# Idents(myeloid.oligo) <- "Miller_clusters"
# p1 = DotPlot(myeloid.astro, features = dot.genes.mac, group.by = "Miller_clusters", 
#              col.min = -0.5, col.max = 1.5, dot.scale = 8, idents = c(0,1,2,3,4)) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100))
# p2 = DotPlot(myeloid.oligo, features = dot.genes.mac, group.by = "Miller_clusters", 
#              col.min = -0.5, col.max = 1.5, dot.scale = 8, idents = c(0,1,2,3)) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100))
# p1 + p2 
# ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_by_cancer.pdf"), width = 17, height = 3.5)
# 
# 
# p1 = DotPlot(myeloid.astro, features = dot.genes.mac, group.by = "grade", 
#              col.min = -0.5, col.max = 1.5, dot.scale = 8) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100))
# p2 = DotPlot(myeloid.oligo, features = dot.genes.mac, group.by = "grade", 
#              col.min = -0.5, col.max = 1.5, dot.scale = 8) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100))
# p1 + p2 
# ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_by_cancer_by_grade.pdf"), width = 17, height = 3.5)
# 
# 
# 
# #for soroban
# cancer = "oligo"
# data.dir = paste0("/mnt/beegfs/home/alyne/stagein/scenicplus_",cancer)
# setwd(data.dir)
# 
# 
# library(igraph)
# library(RColorBrewer)
# library(gage)
# library(annotate)
# library(clusterProfiler)
# library(fgsea)
# library(grid)
# library(readxl)
# 
# if(cancer=="astro"){
#   cor.thresh = 0.45
# }else{
#   cor.thresh = 0.5
# }
# res = 0.18
# 
# fig.dir = paste0("/mnt/beegfs/home/alyne/stageout/figures_sICA_",cancer,"_",cor.thresh)
# dir.create(fig.dir)
# 
# #ggplot palette
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# 
# 
# #set to "malig" if want to run on malignant cells
# fig.str = "malig"
# oligo = readRDS(paste0(data.dir,"/oligo.ss_integrated_RPCA_",fig.str,"_",cancer,".rds"))
# 
# 
# #read in annotated sICA components
# df.annot = as.data.frame(read_xlsx(paste0(fig.dir,"/top_gene_med_",cancer,"_annot.xlsx")))
# colnames(df.annot) = df.annot[1,]
# df.annot = df.annot[-1,]
# prog.orig.order.oligo = colnames(df.annot)
# 
# 
# 
# 
# 
# cancer = "astro"
# data.dir = paste0("/mnt/beegfs/home/alyne/stagein/scenicplus_",cancer)
# if(cancer=="astro"){
#   cor.thresh = 0.45
# }else{
#   cor.thresh = 0.5
# }
# res = 0.18
# 
# fig.dir = paste0("/mnt/beegfs/home/alyne/stageout/figures_sICA_",cancer,"_",cor.thresh)
# astro = readRDS(paste0(data.dir,"/oligo.ss_integrated_RPCA_",fig.str,"_",cancer,"_doub_rem.rds"))
# 
# #read in annotated sICA components
# df.annot = as.data.frame(read_xlsx(paste0(fig.dir,"/top_gene_med_",cancer,"_annot.xlsx")))
# colnames(df.annot) = df.annot[1,]
# df.annot = df.annot[-1,]
# prog.orig.order.astro = colnames(df.annot)
# 
# #find top scores for each cell
# oligo$top_ICA_oligo = sapply(seq(length(colnames(oligo))), function(x) {prog.orig.order.oligo[which.max(oligo@assays$ICA@layers$data[,x])]})
# astro$top_ICA_astro = sapply(seq(length(colnames(astro))), function(x) {prog.orig.order.astro[which.max(astro@assays$ICA@layers$data[,x])]})
# 
# oligo$ICA_mes = oligo$ICA_score_states
# oligo$ICA_mes[which(oligo$top_ICA_oligo=="Reactive astrocyte")] = "Reactive"
# 
# astro$ICA_mes = astro$ICA_score_states
# astro$ICA_mes[which(astro$top_ICA_astro=="Reactive astrocyte")] = "Reactive"
# 
# oligo$ICA_mes_grade = paste("grade",oligo$grade,oligo$ICA_mes)
# astro$ICA_mes_grade = paste("grade",astro$grade,astro$ICA_mes)
# 
# astro$ICA_mes_grade = factor(astro$ICA_mes_grade, levels = c(c("grade 2 astro-like","grade 3 astro-like","grade 4 astro-like",
#                                                                "grade 2 Reactive","grade 3 Reactive","grade 4 Reactive"),
#                                                              unique(astro$ICA_mes_grade)[!(unique(astro$ICA_mes_grade) %in% c("grade 2 astro-like","grade 3 astro-like","grade 4 astro-like",
#                                                                                                                               "grade 2 Reactive","grade 3 Reactive","grade 4 Reactive"))]))
# oligo$ICA_mes_grade = factor(oligo$ICA_mes_grade, levels = c(c("grade 2 astro-like","grade 3 astro-like",
#                                                                "grade 2 Reactive","grade 3 Reactive"),
#                                                              unique(oligo$ICA_mes_grade)[!(unique(oligo$ICA_mes_grade) %in% c("grade 2 astro-like","grade 3 astro-like",
#                                                                                                                               "grade 2 Reactive","grade 3 Reactive"))]))
# 
# dot.genes.astro = c("CD44","TNC","FN1","LAMA4","LTBP1","PTPRM","FZD3","FZD5","LRP5","LRP6")
# 
# 
# fig.dir = '/mnt/beegfs/home/alyne/stageout/figures_cohort_summary'
# Idents(astro) <- "ICA_mes_grade"
# Idents(oligo) <- "ICA_mes_grade"
# p1 = DotPlot(astro, features = dot.genes.astro, group.by = "ICA_mes_grade", idents = c("grade 2 astro-like","grade 3 astro-like","grade 4 astro-like",
#                                                                                        "grade 2 Reactive","grade 3 Reactive","grade 4 Reactive"),
#              col.min = -0.5, col.max = 1.5, dot.scale = 8) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
#   scale_colour_viridis_c(option = "magma", direction=-1)
# p2 = DotPlot(oligo, features = dot.genes.astro, group.by = "ICA_mes_grade", idents = c("grade 2 astro-like","grade 3 astro-like",
#                                                                                        "grade 2 Reactive","grade 3 Reactive"),
#              col.min = -0.5, col.max = 1.5, dot.scale = 8) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
#   scale_colour_viridis_c(option = "magma", direction=-1)
# p1 + p2 #+ plot_layout(guides = "collect")
# ggsave(paste0(fig.dir,"/dotplot_all_malig_cells_by_cancer_by_grade.pdf"), width = 15, height = 3.5)
# 
# 
# 
# Idents(astro) <- "ICA_mes"
# Idents(oligo) <- "ICA_mes"
# p1 = DotPlot(astro, features = dot.genes.astro, group.by = "ICA_mes", idents = c("Reactive","astro-like"),
#              col.min = -0.5, col.max = 1.5, dot.scale = 8) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
#   scale_colour_viridis_c(option = "magma", direction=-1)
# p2 = DotPlot(oligo, features = dot.genes.astro, group.by = "ICA_mes", idents = c("Reactive","astro-like"),
#              col.min = -0.5, col.max = 1.5, dot.scale = 8) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
#   scale_colour_viridis_c(option = "magma", direction=-1)
# p1 + p2 #+ plot_layout(guides = "collect")
# ggsave(paste0(fig.dir,"/dotplot_all_malig_cells_by_cancer.pdf"), width = 15, height = 3.5)
# 
# 
# 
# 
# dot.genes.astro = c("CD44","TNC","FN1","LAMA4","LTBP1","PTPRM","FZD3","FZD5","LRP5","LRP6","PDGFRA","KLRC1","KLRC2","RIMS2","ANGPTL2")
# 
# 
# p1 = DotPlot(astro, features = dot.genes.astro, group.by = "ICA_score_states",
#              col.min = -0.5, col.max = 1.5, dot.scale = 8) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
#   scale_colour_viridis_c(option = "magma", direction=-1)
# p2 = DotPlot(oligo, features = dot.genes.astro, group.by = "ICA_score_states",
#              col.min = -0.5, col.max = 1.5, dot.scale = 8) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
#   scale_colour_viridis_c(option = "magma", direction=-1)
# p1 + p2 #+ plot_layout(guides = "collect")
# ggsave(paste0(fig.dir,"/dotplot_all_malig_cells_by_cancer_by_ICA_score_state.pdf"), width = 18, height = 3.5)
# 
# 
# 
# oligo$ICA_score_state_grade = paste("grade",oligo$grade,oligo$ICA_score_states)
# astro$ICA_score_state_grade = paste("grade",astro$grade,astro$ICA_score_states)
# Idents(astro) <- "ICA_score_state_grade"
# Idents(oligo) <- "ICA_score_state_grade"
# 
# p1 = DotPlot(astro, features = dot.genes.astro, group.by = "ICA_score_state_grade",
#              col.min = -0.5, col.max = 1.5, dot.scale = 8) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
#   scale_colour_viridis_c(option = "magma", direction=-1)
# p2 = DotPlot(oligo, features = dot.genes.astro, group.by = "ICA_score_state_grade",
#              col.min = -0.5, col.max = 1.5, dot.scale = 8) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
#   scale_colour_viridis_c(option = "magma", direction=-1)
# p1 + p2 #+ plot_layout(guides = "collect")
# ggsave(paste0(fig.dir,"/dotplot_all_malig_cells_by_cancer_by_ICA_score_state_by_grade.pdf"), width = 18, height = 3.5)
# 
# 
# wnt.genes = c("LRP5","LRP6","FZD1","FZD3","FZD7",
#               "CTNNB1","AXIN2","LEF1","TCF7",
#               "ROR1","ROR2","RYK","FZD2","FZD5",
#               "VANGL1","VANGL2","PRICKLE1")
# 
# p1 = DotPlot(astro, features = wnt.genes, group.by = "ICA_score_states",
#              col.min = -0.5, col.max = 1.5, dot.scale = 8) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
#   scale_colour_viridis_c(option = "magma", direction=-1)
# p2 = DotPlot(oligo, features = wnt.genes, group.by = "ICA_score_states",
#              col.min = -0.5, col.max = 1.5, dot.scale = 8) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
#   scale_colour_viridis_c(option = "magma", direction=-1)
# p1 + p2 #+ plot_layout(guides = "collect")
# ggsave(paste0(fig.dir,"/dotplot_wnt_genes_all_malig_cells_by_cancer_by_ICA_score_state.pdf"), width = 20, height = 3.5)
# 
# 
# oligo$ICA_score_state_grade = paste("grade",oligo$grade,oligo$ICA_score_states)
# astro$ICA_score_state_grade = paste("grade",astro$grade,astro$ICA_score_states)
# Idents(astro) <- "ICA_score_state_grade"
# Idents(oligo) <- "ICA_score_state_grade"
# 
# 
# p1 = DotPlot(astro, features = wnt.genes, group.by = "ICA_score_state_grade",
#              col.min = -0.5, col.max = 1.5, dot.scale = 8) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
#   scale_colour_viridis_c(option = "magma", direction=-1)
# p2 = DotPlot(oligo, features = wnt.genes, group.by = "ICA_score_state_grade",
#              col.min = -0.5, col.max = 1.5, dot.scale = 8) + 
#   RotatedAxis() + theme_bw() + scale_size(range = c(0, 6), limits = c(0, 100)) +
#   scale_colour_viridis_c(option = "magma", direction=-1)
# p1 + p2 #+ plot_layout(guides = "collect")
# ggsave(paste0(fig.dir,"/dotplot_wnt_genes_all_malig_cells_by_cancer_by_ICA_score_state_by_grade.pdf"), width = 22, height = 3.5)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #make new assay with relevant scores
# my.scores = all.myeloid@meta.data[all.myeloid$disease=="astro",c("Microglia","Complement_Immunosuppressive","Neutrophil","Scavenger_Immunosuppressive","G2_M","IFN_Response",
#                                                                  "Inflammatory_microglia","Macrophage","IL1B_CD83_TNF_Inflammatory","cDC","G1_S","Hypoxia","HS_UPR","Monocyte")]
# 
# assay.v5 <- CreateAssay5Object(data = scale(t(my.scores)))
# all.myeloid_slim <- CreateSeuratObject(assay.v5)
# 
# all.genes = rownames(all.myeloid_slim)
# 
# #also dim reduction
# colnames(my.scores) = paste0("SIXD_",1:ncol(my.scores))
# all.myeloid_slim[["sixD"]] <- CreateDimReducObject(
#   embeddings = as.matrix(my.scores),
#   key = "SIXD_",
#   assay = "RNA"
# )
# 
# 
# all.myeloid_slim <- FindNeighbors(all.myeloid_slim, reduction = "sixD",
#                                   dims = 1:ncol(my.scores), graph.name = c("ICA_state_graph_nn","ICA_state_graph_snn")
# )
# 
# #cluster
# res = 0.45
# all.myeloid_slim <- FindClusters(all.myeloid_slim, resolution = res, cluster.name = paste0("ICA_state_clusters",res), 
#                                  graph = "ICA_state_graph_nn", group.singletons = TRUE)
# #umap
# all.myeloid_slim = RunUMAP(all.myeloid_slim, features = all.genes,
#                            reduction.name = "UMAP_ICA_state", assay="RNA", slot="data")
# 
# FeaturePlot(all.myeloid_slim, features=all.genes, reduction = "UMAP_ICA_state")
# ggsave(paste0(fig.dir,"/umaps_myeloid_cells_all_Miller_astro.pdf"),width=25, height=25)
# 
# 
# #add annotations to plot
# all.myeloid_slim$disease = all.myeloid$disease[which(all.myeloid$disease=="astro")]
# all.myeloid_slim$myeloid_CellType = all.myeloid$myeloid_CellType[which(all.myeloid$disease=="astro")]
# all.myeloid_slim$label_program = all.myeloid$label_program[which(all.myeloid$disease=="astro")]
# all.myeloid_slim$grade = all.myeloid$grade[which(all.myeloid$disease=="astro")]
# all.myeloid_slim$orig.ident = all.myeloid$orig.ident[which(all.myeloid$disease=="astro")]
# 
# 
# DimPlot(all.myeloid_slim, group.by=c(paste0("ICA_state_clusters",res),"disease","myeloid_CellType","label_program","grade","orig.ident"), 
#         reduction = "UMAP_ICA_state", shuffle = T)
# ggsave(paste0(fig.dir,"/umaps_myeloid_cells_clusters_astro.pdf"),
#        width=25, height=12)
# 
# all.myeloid.astro = subset(all.myeloid, subset = disease=="astro")
# 
# all.myeloid.astro$Miller_clusters = all.myeloid_slim$ICA_state_clusters0.45
# # all.myeloid.astro = PrepSCTFindMarkers(all.myeloid.astro)
# 
# DefaultAssay(all.myeloid.astro) <- "RNA"
# all.myeloid.astro <- NormalizeData(all.myeloid.astro)
# all.myeloid.astro <- FindVariableFeatures(all.myeloid.astro)
# all.myeloid.astro <- JoinLayers(all.myeloid.astro)
# markers <- FindAllMarkers(all.myeloid.astro, assay = "RNA", group.by = "Miller_clusters")
# 
# 
# markers = FindAllMarkers(all.myeloid.astro, group.by = "Miller_clusters")
# 
# lapply(unique(all.myeloid.astro$Miller_clusters), function(x) {
#  head(markers[markers$cluster==x & markers$avg_log2FC>0,],30)
# })
# 
# #for astro
# # 0 - intermediate between homeostatic and inflammatory
# # 1 - mature immunoregulatory TAMs/or perivascular TAMs (mostly astro_P0027_ini_S0052)
# # 2 - microglia
# # 3 - inflammatory microglia
# # 4 - phagocytic/scavenger TAM (mostly astro_P0026_ini_S0050)
# # 5 - hypoxia (mostly astro_P0026_ini_S0050)
# # 6 - immunoregulatory (mostly astro_P0028_ini_S0054)
# all.myeloid_slim$cluster_names = sapply(all.myeloid_slim$ICA_state_clusters0.45, function(x) {
#   c("intermediate mic","mature TAMs","microglia","inflammatory","scavenger TAMs","hypoxic","immunoregulatory")[as.numeric(x)]
# })
# DimPlot(all.myeloid_slim, group.by=c(paste0("ICA_state_clusters",res),"cluster_names","myeloid_CellType","label_program","grade","orig.ident"), 
#         reduction = "UMAP_ICA_state", shuffle = T)
# ggsave(paste0(fig.dir,"/umaps_myeloid_cells_clusters_astro.pdf"),
#        width=25, height=12)
# 
# 
# 
# 
# 
# 
# 
# 
# #make new assay with relevant scores
# my.scores = all.myeloid@meta.data[all.myeloid$disease=="oligo",c("Microglia","Complement_Immunosuppressive","Neutrophil","Scavenger_Immunosuppressive","G2_M","IFN_Response",
#                                                                  "Inflammatory_microglia","Macrophage","IL1B_CD83_TNF_Inflammatory","cDC","G1_S","Hypoxia","HS_UPR","Monocyte")]
# 
# assay.v5 <- CreateAssay5Object(data = scale(t(my.scores)))
# all.myeloid_slim <- CreateSeuratObject(assay.v5)
# 
# all.genes = rownames(all.myeloid_slim)
# 
# #also dim reduction
# colnames(my.scores) = paste0("SIXD_",1:ncol(my.scores))
# all.myeloid_slim[["sixD"]] <- CreateDimReducObject(
#   embeddings = as.matrix(my.scores),
#   key = "SIXD_",
#   assay = "RNA"
# )
# 
# 
# all.myeloid_slim <- FindNeighbors(all.myeloid_slim, reduction = "sixD",
#                                   dims = 1:ncol(my.scores), graph.name = c("ICA_state_graph_nn","ICA_state_graph_snn")
# )
# 
# #cluster
# res = 0.4
# all.myeloid_slim <- FindClusters(all.myeloid_slim, resolution = res, cluster.name = paste0("ICA_state_clusters",res), 
#                                  graph = "ICA_state_graph_nn", group.singletons = TRUE)
# #umap
# all.myeloid_slim = RunUMAP(all.myeloid_slim, features = all.genes,
#                            reduction.name = "UMAP_ICA_state", assay="RNA", slot="data")
# 
# FeaturePlot(all.myeloid_slim, features=all.genes, reduction = "UMAP_ICA_state")
# ggsave(paste0(fig.dir,"/umaps_myeloid_cells_all_Miller_oligo.pdf"),width=25, height=25)
# 
# 
# #add annotations to plot
# all.myeloid_slim$disease = all.myeloid$disease[which(all.myeloid$disease=="oligo")]
# all.myeloid_slim$myeloid_CellType = all.myeloid$myeloid_CellType[which(all.myeloid$disease=="oligo")]
# all.myeloid_slim$label_program = all.myeloid$label_program[which(all.myeloid$disease=="oligo")]
# all.myeloid_slim$grade = all.myeloid$grade[which(all.myeloid$disease=="oligo")]
# all.myeloid_slim$orig.ident = all.myeloid$orig.ident[which(all.myeloid$disease=="oligo")]
# 
# 
# DimPlot(all.myeloid_slim, group.by=c(paste0("ICA_state_clusters",res),"disease","myeloid_CellType","label_program","grade","orig.ident"), 
#         reduction = "UMAP_ICA_state", shuffle = T)
# ggsave(paste0(fig.dir,"/umaps_myeloid_cells_clusters_oligo.pdf"),
#        width=25, height=12)
# 
# all.myeloid.oligo = subset(all.myeloid, subset = disease=="oligo")
# 
# all.myeloid.oligo$Miller_clusters = all.myeloid_slim$ICA_state_clusters0.4
# # all.myeloid.oligo = PrepSCTFindMarkers(all.myeloid.oligo)
# 
# DefaultAssay(all.myeloid.oligo) <- "RNA"
# all.myeloid.oligo <- NormalizeData(all.myeloid.oligo)
# all.myeloid.oligo <- FindVariableFeatures(all.myeloid.oligo)
# all.myeloid.oligo <- JoinLayers(all.myeloid.oligo)
# markers <- FindAllMarkers(all.myeloid.oligo, assay = "RNA", group.by = "Miller_clusters")
# # markers = FindAllMarkers(all.myeloid.oligo, group.by = "Miller_clusters")
# 
# lapply(unique(all.myeloid.oligo$Miller_clusters), function(x) {
#   head(markers[markers$cluster==x & markers$avg_log2FC>0,],30)
# })
# 
# #for oligo
# # 0 - microglia
# # 1 - inflammatory
# # 2 - inteferon responsive (mostly sample oligo_P0004_ini_S0009)
# # 3 - mature immunoregulatory TAMs
# 
# 
# all.myeloid_slim$cluster_names = sapply(all.myeloid_slim$ICA_state_clusters0.4, function(x) {
#   c("microglia","inflammatory","interferon responsive","immunoregulatory TAM")[as.numeric(x)]
# })
# DimPlot(all.myeloid_slim, group.by=c(paste0("ICA_state_clusters",res),"cluster_names","myeloid_CellType","label_program","grade","orig.ident"), 
#         reduction = "UMAP_ICA_state", shuffle = T)
# ggsave(paste0(fig.dir,"/umaps_myeloid_cells_clusters_oligo.pdf"), width=25, height=12)
