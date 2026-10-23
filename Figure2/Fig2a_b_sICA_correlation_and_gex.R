##script to project cells into ICA space

rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
data.dir = as.character(args[1])
cancer = as.character(args[2])

#for soroban
cancer = "astro"
data.dir = paste0("/mnt/beegfs/home/alyne/stagein/scenicplus_",cancer)

setwd(data.dir)

library(data.table)
library(gplots)
library(ggplot2)
library(dplyr)
library(Seurat)
library(Signac)
library(igraph)
library(RColorBrewer)
library(gage)
library(annotate)
library(clusterProfiler)
library(fgsea)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(readxl)


#define thresholds depending on subtype
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


#load seurat object of just malignant cells for all samples within subtype
#set to "malig" if want to run on malignant cells
fig.str = "malig"

if(cancer=="astro"){
  oligo.ss = readRDS(paste0(data.dir,"/oligo.ss_integrated_RPCA_",fig.str,"_",cancer,"_doub_rem.rds"))
}else{
  oligo.ss = readRDS(paste0(data.dir,"/oligo.ss_integrated_RPCA_",fig.str,"_",cancer,".rds"))
}

oligo = SplitObject(oligo.ss, split.by = "orig.ident")
samps = sapply(oligo, function(x) unlist(strsplit(names(x$orig.ident)[1],"-"))[3])


#read in meta samples (scores for each component for each cell)
metasamples = list()
for(i in 1:length(samps)){
  metasamples[[i]] = fread(paste0("/mnt/beegfs/home/alyne/stagein/sICA_input/metasamples_",samps[i],".tsv.gz"),
                           sep = "\t", header = T, data.table = F)
  rownames(metasamples[[i]]) = metasamples[[i]][,1]
  metasamples[[i]] = metasamples[[i]][,2:ncol(metasamples[[i]])]
  
  cnames = paste0("mg_",seq(0,(ncol(metasamples[[i]])-1)),
                  "_",samps[i])
  colnames(metasamples[[i]]) = sapply(cnames, function(y) {
    paste(unlist(strsplit(y,"_",fixed = T))[c(1,2,6)],collapse="_")
  })
}

# read in meta genes (scores for genes for each component)
scores = list()
for(i in 1:length(samps)){
  scores[[i]] = fread(paste0("/mnt/beegfs/home/alyne/stagein/sICA_input/metagenes_",samps[i],".tsv.gz"),
                      data.table = F, header = T, colClasses = "numeric", sep = "\t")
  rownames(scores[[i]]) = scores[[i]][,1]
  scores[[i]] = scores[[i]][,2:ncol(scores[[i]])]
}

#gene intersect
genes.list = lapply(scores, function(x) {rownames(x)})
#get intersection of all genes
genes.intersect = Reduce(intersect, genes.list)
#extract only genes in intersect
scores.intersect = lapply(scores, function(x) {
  x = x[sapply(rownames(x), function(y) {y %in% genes.intersect}),]
  #reorder x so that rows are same order as gene list
  x[match(genes.intersect,rownames(x)),]
})

#correlate geps to find which are observed across multiple samples
scores.intersect = lapply(seq(length(scores.intersect)), function(x) {
  cnames = paste0("mg_",seq(0,(ncol(scores.intersect[[x]])-1)),"_",samps[x])
  colnames(scores.intersect[[x]]) = sapply(cnames, function(y) {
    paste(unlist(strsplit(y,"_",fixed = T))[c(1,2,6)],collapse="_")
  })
  return(scores.intersect[[x]])
})
scores.int.all = do.call(cbind,scores.intersect)


#compute correlation
corr.meth = "pearson"
#use all genes
corr = cor(scores.int.all, method = corr.meth)


#read in group membership
df.mem = read.table(paste0(fig.dir,"/sub_graph_members_split_",cancer,".tsv"), sep = "\t", header = T)
df.annot = as.data.frame(read_xlsx(paste0(fig.dir,"/top_gene_med_",cancer,"_annot.xlsx")))
colnames(df.annot) = df.annot[1,]
df.annot = df.annot[-1,]


#give meta components biological names
#re-order programs, lineage programs at top, other programs below
if(cancer=="astro"){
  prog.order = c("Reactive/MES-like","Astro−like1", "OPC−like1",
                 "OPC−like2","OPC−like3","Gliogenesis","Neuro−like",
                 "Cell cycle","Ribosome proteins","Signalling","Stress","Chromatin remodelling",
                 "lnRNA","RNA processing","chr2","chr3","chr4","chr5")
}else{
  prog.order = c("Reactive/MES-like","Astro−like1","Astro−like2","RG/astro synapse","OPC/neuro-like","OPC−like1",
                 "OPC−like2","OPC−like3","Synapse/neuron projection","Neural/developmental",
                 "Cell cycle","Ribosome proteins","MYC targets","Splicing","RNA/proteostatic regulation",
                 "lnRNA","RNA processing","chr1","chr2","chr3","chr4","chr5","chr7")
}


#sort problem with dashes
prog.order <- gsub("\u2212", "-", prog.order)
colnames(df.annot) <- gsub("\u2212", "-", colnames(df.annot))

#get indices to reorder programs
prog.orig.order = colnames(df.annot)
reord.inds = sapply(prog.order, function(x) which(colnames(df.annot)==x))
df.mem = df.mem[,reord.inds]
df.annot = df.annot[,reord.inds]


#filter to take only components in subgroups
subgraphs = colnames(df.mem)
col.pal = gg_color_hue(ncol(df.mem))
plot_corr = corr
cols.corr = sapply(rownames(plot_corr), function(x) {
  #find if component is in one of the subgraphs
  logi = apply(df.mem, 2, function(y) {x %in% y})
  if(any(logi)){
    sg = colnames(df.mem)[which(logi)]
    col.pal[which(subgraphs==sg)]
  }else{
    "white"
  }
})
#remove anything white
plot_corr = plot_corr[which(cols.corr!="white"),which(cols.corr!="white")]
cols.corr = cols.corr[which(cols.corr!="white")]


##############
## make plot of correlation of ICA components (Fig 1A)
## onyl components which end up in metacomponents are included
##############

# # Hierarchical clustering
hc <- hclust(dist(plot_corr), method = "complete")

#assign groups to each row
groups = sapply(rownames(plot_corr), function(x) {
  #find if component is in one of the subgraphs
  logi = apply(df.mem, 2, function(y) {x %in% y})
  colnames(df.annot)[which(logi)]
})

# get order within dendrogram
ord <- hc$order
groups_ord <- groups[ord]
# Make factor with correct order
groups_ord <- factor(groups_ord, levels = prog.order)
#make colours fro each set
group_cols <- setNames(col.pal, prog.order)

col_fun <- colorRamp2(c(-0.7, 0, 0.7), c("indianred1","white","royalblue1"))
row_anno <- rowAnnotation(
  Cluster = anno_block(
    gp = gpar(fill = group_cols[levels(groups_ord)], col = NA),  # col = NA removes border
    labels = paste(levels(groups_ord)),
    labels_gp = gpar(fontsize = 6),
    labels_rot = 0
  ),
  width = unit(3.5, "cm"),
  show_annotation_name = FALSE
)
ht <- Heatmap(
  plot_corr[ord,ord],
  name = "Correlation",
  col = col_fun,
  cluster_rows = F,
  cluster_columns = F,
  row_split = groups_ord,
  column_split = groups_ord,
  row_gap = unit(0, "mm"),
  column_gap = unit(0, "mm"),
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  left_annotation = row_anno,
  border=F,
  show_row_names = F,
  show_column_names = F,
  row_title = NULL,
  column_title = NULL,
  rect_gp = gpar(col = NA),
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(
      x, y,
      w + unit(0.3, "pt"),
      h + unit(0.3, "pt"),
      gp = gpar(fill = fill, col = NA)
    )
  }
)
pdf(paste0(fig.dir,"/final_heatmap_pearson_correlation_filtered_",cancer,".pdf"), width = 17, height = 14, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
decorate_heatmap_body("Correlation", {
  n <- length(unique(groups_ord))
  for(i in seq_len(n)) {
    seekViewport(paste0("Correlation_heatmap_body_", i, "_", i))
    grid.rect(x = 0.5, y = 0.5,width = 1, height = 1,gp = gpar(col = "black", lty = 2, lwd = 2, fill = NA))
    upViewport()
  }
})
dev.off()


################
## Correlation of scores across cells i.e. which programs tend to be mutually expressed
## Fig 1E
################

to_corr = as.matrix(t(oligo.ss@assays$ICA@layers$data[,]))
colnames(to_corr) = prog.orig.order
M = cor(to_corr, method="pearson")
col_fun <- colorRamp2(c(-0.55, 0, 0.55), c("blue","white","red"))
ht = Heatmap(as.matrix(M), name = "Correlation",
             col = col_fun, row_split = 3, column_split = 3,
             clustering_distance_rows = "pearson",
             clustering_distance_columns = "pearson",
             clustering_method_rows = "ward.D2",
             clustering_method_columns = "ward.D2")
pdf(paste0(fig.dir,"/corr_ICA_components_",cancer,"_pearson.pdf"), width = 8, height = 7, useDingbats = FALSE)
draw(ht)
dev.off()



####################
### plot gene expression heatmap
####################

#heatmap showing expression values for top genes for each component for subsample of cells within subtype (Fig 2B)
#each row is a gene, split into groups according to component
#each column a cell, organised according to top dimension
#for each program, also show proportion of cells

Nplot.gene = 40
Nplot.cell = 10000
cells.to.samp = sample(ncol(oligo.ss@assays$SCT@data),Nplot.cell,replace=F)
unscaled.mat = oligo.ss@assays$SCT@data

#####
#make data frame for column annotations
#####

col_anno_df = data.frame(oligo.ss@meta.data[cells.to.samp,colnames(oligo.ss@meta.data) %in% 
                                              c("ICA_score_states","ICA_cell_states","cyc.status","atlas.label2","orig.ident")])
#make list of colours for each category
cols.list.col = lapply(seq(ncol(col_anno_df)), function(x) {
  cols = gg_color_hue(length(unique(col_anno_df[,x])))
  names(cols) = sort(unique(col_anno_df[,x]))
  return(cols)
})
names(cols.list.col) = colnames(col_anno_df)

#add which ICA component has highest projected score for each cell
col_anno_df$top_ICA = sapply(cells.to.samp, function(x) {prog.orig.order[which.max(oligo.ss@assays$ICA@layers$data[,x])]})
col_anno_df$top_ICA = factor(col_anno_df$top_ICA, levels = prog.order)
cols.list.col$top_ICA = group_cols

#also find top module score for each cell
to.plot.sc = t(oligo.ss@meta.data[cells.to.samp,which(grepl("^ICA_.+_1$",colnames(oligo.ss@meta.data)))])
rownames(to.plot.sc) = prog.orig.order

#add column to annotation df of top module score
col_anno_df$top_ICA_mod = sapply(rownames(col_anno_df), function(x) {
  rownames(to.plot.sc)[which.max(to.plot.sc[,which(colnames(to.plot.sc)==x)])]
})
col_anno_df$top_ICA_mod = factor(col_anno_df$top_ICA_mod, levels = prog.order)
cols.list.col$top_ICA_mod = cols.list.col$top_ICA

#add patient and sample status (ini/rec)
col_anno_df$patient = sapply(rownames(col_anno_df), function(x) {unlist(strsplit(x,"_"))[2]})
col_anno_df$ini_rec = sapply(rownames(col_anno_df), function(x) {unlist(strsplit(x,"_"))[3]})


#####
#get genes to plot
#####

genes.to.plot = lapply(prog.order, function(x) {df.annot[1:Nplot.gene,which(colnames(df.annot)==x)]})
gene.inds = lapply(genes.to.plot, function(x) {
  sapply(x, function(y) {which(rownames(unscaled.mat)==y)})
})

#vector to split heatmap between each metaprogram
split_groups = unlist(lapply(prog.order, function(x) {rep(x,Nplot.gene)}))
split_groups = factor(split_groups, levels = prog.order)

#get values from SCT data matrix
vals.to.plot = lapply(seq(length(gene.inds)), function(x) {
  rel.rows = unscaled.mat[gene.inds[[x]],cells.to.samp]
  #add suffix so can plot same genes multiple times
  rownames(rel.rows) = paste0(rownames(rel.rows),".",x)
  return(rel.rows)
})
data.mat = do.call("rbind",vals.to.plot)

#scale
to.plot <- t(scale(t(data.mat)))
to.plot[is.na(to.plot)] <- 0



# genes to label for each metaprogram
if(cancer=="astro"){
  sel.list = list(Reactive_astrocyte = c("GFAP","TNC","CD44"),
                  Astro_like2 = c("SLC1A3","SPARCL1","NFIA"),
                  OPC_like1 = c("VCAN","SOX6","LHFPL3","FGF14"),
                  OPC_like2 = c("PDGFRA","SOX5","CNTN1"),
                  OPC_like3 = c("EGFR","NRXN1","SLIT1"),
                  NeuroGlio_genesis = c("ZEB2","JAG1","NPAS3"),
                  Neuro_like = c("DCC","MYT1L","ROBO1","ROBO2"),
                  Cell_cycle = c("TOP2A", "CENPF"),
                  Ribosome_proteins = c("EEF1A1", "RPLP0", "RPL13A"),
                  Signalling = c("VAV3","TIAM1","CHRM3"),
                  Stress = c("FOS", "JUN", "JUNB", "EGR1"),
                  Chromatin_remodelling = c("KDM4B","GATAD2B","BPTF","SMARCA4"),
                  lnRNA = c("MALAT1","AQP4-AS1","MIR100HG"),
                  mRNA_processing = c("CNOT1","FIP1L1","AQR"),
                  chr2 = character(),chr3 = character(),chr4 = character(),chr5 = character())
}else{
  sel.list = list(Reactive_astrocyte = c("GFAP","CD44","TNC"),
                  Astro_like1 = c("GFAP","APOE","CLU"),
                  Astro_like2 = c("KCNN3","SPARCL1"),
                  RGastro_synapse = c("SLC1A3","SPARCL1","GRIA1","NRXN3"),
                  OPCneuro_like = c("OPCML","SOX5","NRXN1","NLGN1"),
                  OPC_like1 = c("VCAN","SOX6","LHFPL3"),
                  OPC_like2 = c("PDGFRA","PTPRZ1","CNTN1"),
                  OPC_like3 = c("EGFR","SLIT1"),
                  SynapseNeuron_projection = c("GRIA3",'ERBB4',"NRCAM"),
                  Neural_developmental = c("NFIX","GRIK4","DSCAML1"),
                  Cell_cycle = c("EZH2", "NUSAP1"),
                  Ribosome_proteins = c("EEF1A1", "RPLP0", "RPL13A"),
                  MYC_targets = c("MTHFD1L","FARSB"),
                  Splicing = c("PNISR","LUC7L3","RBM25"),
                  RNA_proteostatic_regulation = c("ARID2", "ARID4B", "SETD2"),
                  lnRNA = c("LHFPL3-AS1","AQP4-AS1","MALAT1"),
                  RNA_processing = c("PRPF6", "DDX3X", "DDX42", "SAFB2"),
                  chr1 = character(),chr2 = character(),chr3 = character(),chr4 = character(),chr5 = character(),chr7 = character())
}
sel <- unlist(sel.list)

#find the correct names in the matrix with suffices in the correct program
genes.suff = unlist(lapply(seq(length(sel.list)), function(x) {
  sapply(sel.list[[x]], function(y) {
    rownames(vals.to.plot[[x]])[which(grepl(paste0(y,"."), rownames(vals.to.plot[[x]]),fixed=T))]
  })
}))


# show proportion of cells in each sample which has this program as highest score
#for each cell, find which component has highest score
max.ICA.inds = apply(oligo.ss@assays$ICA@layers$data, 2, which.max)
#for each program, find proportion of cells with highest in each sample
prop.per.comp = lapply(seq(nrow(oligo.ss@assays$ICA@layers$data)), function(x) {
  sapply(unique(oligo.ss$orig.ident), function(y) {
    all.max.inds = max.ICA.inds[which(oligo.ss$orig.ident==y)]
    sum(all.max.inds==x)/length(all.max.inds)
  })
})
prop.per.comp = prop.per.comp[reord.inds]
df.prop.comps = data.frame(component = unlist(lapply(prog.order, function(x) {rep(x,length(unique(oligo.ss$orig.ident)))})),
                           proportion = unlist(prop.per.comp))
df.prop.comps$component = factor(df.prop.comps$component, levels = prog.order)



######################
#heatmap with no clustering and ordered by highest ICA component in projection
#keep only components which have clear populations
######################

if(cancer=="astro"){
  progs.rem = c("chr2","chr3","chr4","chr5")
}else{
  progs.rem = c("chr1","chr2","chr3","chr4","chr5","chr7")
}


prog.order.filt = prog.order[!(prog.order %in% progs.rem)]
col.pal.filt = gg_color_hue(length(prog.order.filt))
group_cols.filt <- setNames(col.pal.filt, prog.order.filt)

to.plot.filt = to.plot[!(split_groups %in% progs.rem),]

######## shouod this be group_cols.filt
row_anno <- rowAnnotation(
  Cluster = anno_block(gp = gpar(fill = group_cols[prog.order.filt]), labels = prog.order.filt,
                       labels_gp = gpar(fontsize = 8), labels_rot = 0),
  width = unit(3.5, "cm"), show_annotation_name = FALSE
)
#remove unrequired components
df.prop.comps.filt = df.prop.comps[!(df.prop.comps$component %in% progs.rem),]
df.prop.comps.filt$component = droplevels(df.prop.comps.filt$component)


#get indices to label
idx.filt <- which(rownames(to.plot.filt) %in% genes.suff)


#et row splits only for retained programs
split_groups.filt = split_groups[!(split_groups %in% progs.rem)]
split_groups.filt = droplevels(split_groups.filt)

#make right annotation showing proportion of cells with that metacomponent as its max projection (boxplot across samples)
prop_list.filt <- split(df.prop.comps.filt$proportion, df.prop.comps.filt$component)
right_anno <- rowAnnotation(
  mark = anno_mark(
    at = which(rownames(to.plot.filt) %in% genes.suff),
    labels = sapply(rownames(to.plot.filt)[idx.filt], function(x) strsplit(x, ".", fixed=TRUE)[[1]][1]),
    labels_gp = gpar(fontsize=7)
  ),
  proportion = anno_block(
    gp = gpar(fill = NA),
    panel_fun = function(index, nm) {
      vals <- prop_list.filt[[nm]]
      vals <- vals[is.finite(vals)]
      if(length(vals) < 2) return(NULL)
      q <- quantile(vals, c(0, .25, .5, .75, 1))
      pushViewport(viewport(xscale = c(0, max(unlist(prop_list.filt))), yscale = c(0,1)))
      grid.rect(x=(q[2]+q[4])/2, y=0.5, width=q[4]-q[2], height=0.5, gp=gpar(fill="grey80", col="grey40"))
      grid.lines(x=c(q[3],q[3]), y=c(0.25,0.75), gp=gpar(lwd=2))
      grid.lines(x=c(q[1],q[2]), y=c(0.5,0.5))
      grid.lines(x=c(q[4],q[5]), y=c(0.5,0.5))
      # Draw x-axis only for the last group (bottom)
      if(nm == tail(names(prop_list.filt),1)) {
        grid.xaxis(at = pretty(c(0, 0.4)), gp=gpar(fontsize=9))
      }
      popViewport()
    },
    width = unit(4.2, "cm")
  ), gap = unit(20, "mm")
)


#order columns by top ICA component
to.plot.filt.reord = to.plot.filt[,order(col_anno_df$top_ICA)]
col_anno_df.reord = col_anno_df[order(col_anno_df$top_ICA),]
top_anno <- HeatmapAnnotation(df = col_anno_df.reord[!(col_anno_df.reord$top_ICA %in% progs.rem),c("cyc.status","top_ICA")], 
                              col = cols.list.col,show_legend = c(TRUE,FALSE))

#finally plot heatmap
col_fun <- colorRamp2(c(-1.6, 0, 1.6), c("royalblue1", "white", "darkred"))
ht <- Heatmap(
  as.matrix(to.plot.filt.reord[,!(col_anno_df.reord$top_ICA %in% progs.rem)]),
  col = col_fun,
  name = "Scaled gene expression",
  cluster_rows = F,
  cluster_columns = F,
  row_split = split_groups.filt,       # vector of cluster IDs
  column_split = col_anno_df.reord$top_ICA[!(col_anno_df.reord$top_ICA %in% progs.rem)],       # vector of cluster IDs
  row_gap = unit(0.3, "mm"),
  column_gap = unit(0.3, "mm"),
  bottom_annotation  = top_anno,
  left_annotation = row_anno,
  right_annotation = right_anno,
  show_row_names = F,
  show_column_names = F,
  row_title = NULL,
  column_title = NULL,
  show_column_dend = F,
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect( x, y, w + unit(0.3, "pt"), h + unit(0.3, "pt"),
      gp = gpar(fill = fill, col = NA) )
  }
)
pdf(paste0(fig.dir,"/final_heatmap_gex_filt_nocluster_",cancer,".pdf"), width = 15, height = 7, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",merge_legends = TRUE)
dev.off()






