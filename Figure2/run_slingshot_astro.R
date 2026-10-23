
library(data.table)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(slingshot)
library(mclust)
library(grDevices)
library(cowplot)
library(scater)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(destiny)
library(scran)
library(purrr)
library(uwot)
library(readxl)
library(viridis)
library(dplyr)
library(forcats)
library(stringr)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(ggrastr)

pal = viridis(n=50)

#set to "_rem" if want to remove samples, otherwise
rem = ""
samps.to.rem = c("astro_P0028_ini_S0054")

dr = "dm_ica"

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols.TF = gg_color_hue(2)

data.dir = "/Users/alyne/Library/CloudStorage/Dropbox/Curie/Curie_StCloud/anne-marie/10x_RNA_ATAC_oligo_analysis/trajectory_analysis_local"
pt.sz = 0.1
cancer = "astro"
j="all"

setwd(paste0(data.dir,"/",cancer,"_ss_rem_cc_",j))
fig.dir = paste0("figures_dr_",cancer,"_",j)
dir.create(fig.dir)

source(paste0("../../corrected_GSEA_function.R"))

#read in SCT data matrix for this subsample
data.sct = fread(paste0(cancer,"_samps_ss_full_data_rem_cc_",j,".tsv.gz"), data.table=F)
rownames(data.sct) = data.sct[,1]
data.sct = data.sct[,2:ncol(data.sct)]

#read in cell metadata
cell.md = fread(paste0("metadata_",cancer,"_samps_ss_rem_cc_",j,".tsv.gz"), data.table=F)
rownames(cell.md) = colnames(data.sct)

#read in ICA metadata
cell.md.ICA = fread(paste0("metadata_ICA_",cancer,"_samps_ss_rem_cc_",j,".tsv.gz"), data.table=F)
rownames(cell.md.ICA) = cell.md.ICA[,1]
cell.md.ICA = cell.md.ICA[,2:ncol(cell.md.ICA)]

#read in sample summary
summ.file = paste0("/Users/alyne/Library/CloudStorage/Dropbox/Curie/Curie_StCloud/anne-marie/data/",
                   "our_data/sample_summaries/sample_summary_astro_oligo_anon.xlsx")
samp.md = as.data.frame(read_excel(summ.file))


# #read in projection into ICA space
data.ica = fread(paste0(cancer,"_samps_ss_rem_cc_",j,".tsv.gz"), data.table=F)
rownames(data.ica) = data.ica[,1]
data.ica = data.ica[,2:ncol(data.ica)]

if(rem=="_rem"){
  if(length(samps.to.rem)>0){
    for(i in 1:length(samps.to.rem)){
      data.sct = data.sct[,which(!grepl(samps.to.rem[i],colnames(data.sct)))]
      data.ica = data.ica[,which(!grepl(samps.to.rem[i],colnames(data.ica)))]
      cell.md = cell.md[which(!grepl(samps.to.rem[i],rownames(cell.md))),]
      cell.md.ICA = cell.md.ICA[which(!grepl(samps.to.rem[i],rownames(cell.md.ICA))),]
    }
  }
}


#make sce object
#actually these are more like log counts from data slot of SCT assay - need to check this is ok for fitGAM
sce <- SingleCellExperiment(assays = list(counts = as(data.sct, "sparseMatrix")),
                            colData = cell.md)


#compute diffusion map
set.seed(12345)
dm <- DiffusionMap(t(data.ica),n_pcs=NULL)
reducedDim(sce, "DM") <- cbind(dm$DC1,dm$DC2)



df = data.frame(DR1 = dm$DC1, DR2 = dm$DC2, DR3 = dm$DC3)

#k means clustering
if(cancer=="astro"){
  ncent = 8
}else if(cancer=="oligo"){
  ncent = 7
}

cl2 <- kmeans(df, centers = ncent)$cluster
colData(sce)$kmeans <- cl2
#kmeans and metadata
df = cbind(df, cell.md, kmeans_DR = factor(colData(sce)$kmeans))

#same for umap projection
cl2 <- kmeans(df[,which(colnames(df) %in% c("UMAP1","UMAP2"))], centers = ncent)$cluster
#kmeans and metadata
df = cbind(df, kmeans_UMAP = factor(cl2))

#df for ica component scores
df.ICA = cbind(data.frame(DR1 = dm$DC1, DR2 = dm$DC2), UMAP1 = df$UMAP1, UMAP2 = df$UMAP2, cell.md.ICA)
df.ICA2 = df.ICA[sample(seq(nrow(df.ICA)),nrow(df.ICA)),]

for(dr in c("DR","UMAP")){
  plot.list = list()
  count = 1
  #plot LLE with clusters
  plot.list[[count]] = ggplot(df, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour=.data[[paste0("kmeans_",dr)]])) + geom_point(size=pt.sz) + 
    guides(color = guide_legend(override.aes = list(size=2)))
  count = count + 1
  #tirosh cell types
  df2 = df[sample(seq(nrow(df)),nrow(df)),]
  plot.list[[count]] = ggplot(df2, aes(x = DR1, y = DR2, colour = ICA_cell_states)) + geom_point(size=pt.sz) + 
    guides(color = guide_legend(override.aes = list(size=2)))
  count = count + 1
  plot.list[[count]] = ggplot(df2, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour = ICA_score_states)) + geom_point(size=pt.sz) + 
    guides(color = guide_legend(override.aes = list(size=2)))
  count = count + 1
  plot.list[[count]] = ggplot(df2, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour = cell.type.sr.imm)) + geom_point(size=pt.sz) + 
    theme(legend.key.size = unit(0.2, 'cm')) + guides(colour = guide_legend(override.aes = list(size=3)))
  count = count + 1
  #SingleR cell type
  #sort df so neurons are at bottom
  df2 = rbind(df[!(df$cell.type.sr.imm %in% c("GABAergic","Glutamatergic")),],
              df[df$cell.type.sr.imm %in% c("GABAergic","Glutamatergic"),])
  plot.list[[count]] = ggplot(df2, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour = cell.type.sr.imm)) + geom_point(size=pt.sz) + 
    theme(legend.key.size = unit(0.2, 'cm')) + guides(colour = guide_legend(override.aes = list(size=3)))
  count = count + 1
  #proportion cycling for each cluster and cycling cells
  prop.cyc.cl = sapply(unique(df[[paste0("kmeans_",dr)]]), function(x) {
    sum(df$cyc.status[which(df[[paste0("kmeans_",dr)]]==x)]=="cycling")/length(df$cyc.status[which(df[[paste0("kmeans_",dr)]]==x)])
  })
  df$prop.cyc = sapply(df[[paste0("kmeans_",dr)]], function(x) {
    prop.cyc.cl[which(unique(df[[paste0("kmeans_",dr)]])==x)]
  })
  plot.list[[count]] = ggplot(df, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour = prop.cyc)) + geom_point(size=pt.sz)
  count = count + 1
  df2 = rbind(df[!(df$cyc.status %in% c("cycling")),], df[df$cyc.status %in% c("cycling"),])
  plot.list[[count]] = ggplot(df2, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour = cyc.status)) + geom_point(size=pt.sz) + 
    guides(color = guide_legend(override.aes = list(size=2)))
  count = count + 1
  #plot cytotrace value
  df2 = df[sample(seq(nrow(df)),nrow(df)),]
  plot.list[[count]] = ggplot(df2, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour = CytoTRACE)) + geom_point(size=pt.sz)
  count = count + 1
  #sample and ini/rec and icnv score
  plot.list[[count]] = ggplot(df2, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour = sample)) + geom_point(size=pt.sz) + 
    theme(legend.key.size = unit(0.2, 'cm')) + guides(colour = guide_legend(override.aes = list(size=3)))
  count = count + 1
  plot.list[[count]] = ggplot(df2, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour = ini_rec)) + geom_point(size=pt.sz) + 
    guides(color = guide_legend(override.aes = list(size=2)))
  count = count + 1
  plot.list[[count]] = ggplot(df2, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour = icnv.score)) + geom_point(size=pt.sz)+ 
    scale_colour_gradientn(colours = pal, na.value = "grey90")
  count = count + 1
  plot.list[[count]] = ggplot(df2, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour = as.factor(grade))) + geom_point(size=pt.sz) + 
    guides(color = guide_legend(override.aes = list(size=2)))
  count = count + 1
  plot.list[[count]] = ggplot(df2, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour = treatment)) + geom_point(size=pt.sz) + 
    guides(color = guide_legend(override.aes = list(size=2)))
  count = count + 1
  plot.list[[count]] = ggplot(df2, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour = atlas.label)) + geom_point(size=pt.sz) + 
    theme(legend.key.size = unit(0.2, 'cm')) + guides(colour = guide_legend(override.aes = list(size=3)))
  count = count + 1
  plot.list[[count]] = ggplot(df2, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour = gest_week)) + geom_point(size=pt.sz) + 
    guides(color = guide_legend(override.aes = list(size=2)))
  count = count + 1
  p1 = wrap_plots(plot.list, ncol=4)
  print(p1)
  ggsave(paste0(fig.dir,"/pre_traj_",cancer,"_",j,"_",dr,rem,".pdf"), width = 20, height = 12)
  
  # #plot version with ICA scores
  # plot.list2 = list()
  # for(i in 1:ncol(cell.md.ICA)){
  #   plot.list2[[i]] = ggplot(df.ICA2, aes(x=.data[[paste0(dr,"1")]], y=.data[[paste0(dr,"2")]], colour = .data[[colnames(cell.md.ICA)[i]]])) + 
  #     geom_point(size=pt.sz) + scale_color_viridis()
  # }
  # p2 = wrap_plots(plot.list2, ncol=4)
  # print(p2)
  # ggsave(paste0(fig.dir,"/pre_traj_ICA_",cancer,"_",j,"_",dr,rem,".pdf"), width = 20, height = 17)
}

save(df, df.ICA, sce, file = paste0("pre_traj_output_",cancer,"_",j,rem,".RData"))




#####################
## Get trajectory with slingshot
#####################

load(paste0("pre_traj_output_",cancer,"_",j,rem,".RData"))

#half cluster 6 so can encourage trajectory to go through main density
sce$kmeans[which(sce$kmeans==5 & df$DR1>0.002)] = 9
#fit trajectory and curves etc
end.cl = c(2,1)
start.cl = 4


df = df[,1:21]
#remove any slingshot 
colData(sce) = colData(sce)[,which(!grepl("^slingPseudotime_",colnames(colData(sce))))]
sce$slingshot = NULL
sce <- slingshot(sce, clusterLabels = 'kmeans', reducedDim = 'DM',
                 start.clus= as.character(start.cl), end.clus = as.character(end.cl))

#add pseudotime to data frame
#count how many pseudotimes there are
Nps = length(slingCurves(sce))
inds = (ncol(colData(sce))-Nps+1):ncol(colData(sce))

for(i in 1:length(inds)){
  df = cbind(df, colData(sce)[,inds[i]])
  colnames(df)[ncol(df)] = paste0("slingPseudotime_",i)
  if(i==1){
    df.crv = data.frame(crv1_1 = slingCurves(sce)[[1]]$s[,1], crv1_2 = slingCurves(sce)[[1]]$s[,2])
  }else{
    df.crv = cbind(df.crv, slingCurves(sce)[[i]]$s[,1], slingCurves(sce)[[i]]$s[,2])
    colnames(df.crv)[ncol(df.crv)-1] = paste0("crv",i,"_1")
    colnames(df.crv)[ncol(df.crv)] = paste0("crv",i,"_2")
  }
}

add_curves <- function(plt, inds){
  for(i in 1:length(inds)){
    plt = plt + geom_path(data = df.crv, aes(x=.data[[paste0("crv",i,"_1")]],
                                             y=.data[[paste0("crv",i,"_2")]]),col="black")
  }
  return(plt)
}

#plot with pseudotime and trajectory
plot.list = list()
plot.list[[1]] = ggplot(df, aes(x = DR1, y = DR2, colour = kmeans_DR)) + geom_point(size=pt.sz) + guides(color = guide_legend(override.aes = list(size=2)))
plot.list[[1]] = add_curves(plot.list[[1]], inds)
plot.list[[2]] = ggplot(df, aes(x = DR1, y = DR2, colour = ICA_score_states)) + geom_point(size=pt.sz) + guides(color = guide_legend(override.aes = list(size=2)))
plot.list[[2]] = add_curves(plot.list[[2]], inds)
for(m in 1:length(inds)){
  plot.list[[2+m]] = ggplot(df, aes(x = DR1, y = DR2, colour = .data[[paste0("slingPseudotime_",m)]])) + geom_point(size=pt.sz)
  plot.list[[2+m]] = add_curves(plot.list[[2+m]], inds)
}
p1 = wrap_plots(plot.list)
ggsave(paste0(fig.dir,"/traj_",cancer,"_",j,"_",dr,"_rem.pdf"), width = 12, height = 8)


save(df, df.ICA, sce, file = paste0("post_traj_output_",cancer,"_",j,rem,".RData"))


load(paste0("post_traj_output_",cancer,"_",j,rem,".RData"))

df = df[,1:23]
df2 = df
Nps = length(slingCurves(sce))
inds = (ncol(colData(sce))-Nps+1):ncol(colData(sce))
#compute genes which correlate well with each trajectory
#and do enrichment
for(i in 1:length(inds)){
  sce2 = sce
  corr_ps = lapply(seq(nrow(sce2)), function(x) {
    res = cor.test(df2[[paste0("slingPseudotime_",i)]],counts(sce2)[x,], 
                   use = "complete.obs", method = "pearson")
    return(c(res$estimate, res$p.value))
  })
  df.corr = do.call("rbind",corr_ps)
  df.corr = as.data.frame(df.corr)
  rownames(df.corr) = rownames(sce2)
  colnames(df.corr) = c("corr_ps","pval_1")
  
  # Nplot = 6
  # for(direc in c("up","down")){
  #   plot.list = list()
  #   plot.list2 = list()
  #   #plot top 5, both expression on trajectory and dot plot of expr v pseudotime
  #   for(b in 1:Nplot){
  #     if(direc=="up"){
  #       gene = rownames(df.corr[order(df.corr$corr_ps,decreasing=T),])[b]
  #     }else{
  #       gene = rownames(df.corr[order(df.corr$corr_ps),])[b]
  #     }
  #     df3 = data.frame(pseudotime = colData(sce2)[,inds[i]], expr = counts(sce2)[which(rownames(sce2)==gene),],
  #                      cluster = df$kmeans_DR)
  #     plot.list2[[b]] = ggplot(df3, aes(x = pseudotime, y = expr, color = cluster)) +
  #       geom_point(size=pt.sz) + labs(title=gene)
  #     
  #     df3 = df2
  #     df3 = cbind(df3, counts(sce2)[which(rownames(sce2)==gene),])
  #     colnames(df3)[length(colnames(df3))] = gene
  #     plot.list[[b]] = ggplot(df3, aes(x = DR1, y = DR2, color = .data[[gene]])) + geom_point(size=pt.sz) +
  #       geom_path(data = df.crv, aes(x=crv1_1,y=crv1_2),col="black") +
  #       geom_path(data = df.crv, aes(x=crv2_1,y=crv2_2),col="black")
  #   }
  #   p1 = wrap_plots(c(plot.list,plot.list2), ncol = Nplot, byrow = T)
  #   print(p1)
  #   ggsave(paste0(fig.dir,"/ps",i,"_gene_expr_",direc,"_",cancer,"_",j,"_DR.pdf"), width = 20, height = 5)
  # }
  
  #enrichment in top
  nUse = 100
  pval.thr = 1
  qval.thr = 1
  for(direc in c("up","down")){
    if(direc=="up"){
      gene.list = rownames(df.corr[order(df.corr$corr_ps,decreasing=T),])[1:nUse]
    }else{
      gene.list = rownames(df.corr[order(df.corr$corr_ps),])[1:nUse]
    }
    eGOBP.high = enrichGO(gene.list, OrgDb ="org.Hs.eg.db", keyType ="SYMBOL", ont ="BP",
                          pvalueCutoff = pval.thr, minGSSize = 10, maxGSSize = 500,
                          qvalueCutoff = qval.thr, universe = rownames(sce2))
    if(sum(eGOBP.high@result$pvalue<pval.thr & eGOBP.high@result$p.adjust<pval.thr &
           eGOBP.high@result$qvalue<qval.thr)>0){
      pdf(paste0(fig.dir,"/ps",i,"_enrichment_GOBP_top_",nUse,
                 "_",cancer,"_",j,"_DR_",direc,".pdf"), height = 3, width = 8)
      print(barplot(eGOBP.high, showCategory = 5, label_format = 35, font.size = 18) +
              theme(legend.text = element_text(size = 14),legend.title = element_text(size = 16)))
      dev.off()
    }
    
    m_t2g = msigdbr(species = "Homo sapiens") %>% dplyr::select(gs_name, gene_symbol)
    eMSIGDB.high = enricher(gene.list, TERM2GENE = m_t2g, universe = rownames(sce2),
                            pvalueCutoff = pval.thr, minGSSize = 10, maxGSSize = 500, qvalueCutoff = qval.thr)
    if(sum(eMSIGDB.high@result$pvalue<pval.thr & eMSIGDB.high@result$p.adjust<pval.thr &
           eMSIGDB.high@result$qvalue<qval.thr)>0){
      pdf(paste0(fig.dir,"/ps",i,"_enrichment_MSIGDB_top_",nUse,
                 "_",cancer,"_",j,"_DR_",direc,".pdf"), height = 4, width = 8)
      print(barplot(eMSIGDB.high, showCategory = 5, label_format = 35, font.size = 18) + 
              theme(legend.text = element_text(size = 14),legend.title = element_text(size = 16)))
      dev.off()
    }
  }
  
  # cats = c("H","C8","C2")
  # gsea.input = df.corr[order(df.corr$corr_ps,decreasing=T),1]
  # names(gsea.input) = rownames(df.corr[order(df.corr$corr_ps,decreasing=T),])
  # barplot.list = list()
  # count = 1
  # for(cat in cats){
  #   if(cat=="C2"){
  #     msigdbr_df <- msigdbr(species = "human", collection = cat, subcollection = "CP:WIKIPATHWAYS")
  #   }else{
  #     msigdbr_df <- msigdbr(species = "human", collection = cat)
  #   }
  #   pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
  #   #remove repeated genes as this was causing a problem in GSEA
  #   pathwaysH = lapply(pathwaysH, unique)
  #   
  #   fgseaRes = fgseaCellIDcorr(pathways=pathwaysH, gsea.input, minSize = 15, maxSize = 500)
  #   
  #   #tidy results
  #   fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
  #   fgseaResPos = fgseaResTidy %>% filter(NES > 0)
  #   fgseaResNeg = fgseaResTidy %>% filter(NES < 0)
  #   
  #   #Plot the normalized enrichment scores.
  #   #Color the bar indicating whether or not the pathway was significant
  #   nPrint = 10
  #   if(nrow(fgseaResPos)>=(nPrint)){
  #     fgseaResPos = fgseaResPos[1:nPrint,]
  #   }
  #   if(nrow(fgseaResNeg)>=(nPrint)){
  #     fgseaResNeg = fgseaResNeg[(nrow(fgseaResNeg)-nPrint+1):nrow(fgseaResNeg),]
  #   }
  #   fgseaResTidy = bind_rows(fgseaResPos, fgseaResNeg)
  #   fgseaResTidy$psign = sign(fgseaResTidy$NES)*(-log10(fgseaResTidy$padj))
  #   fgseaResTidy$significance = fgseaResTidy$padj<0.05
  #   fgseaResTidy$TFcolor = sapply(fgseaResTidy$significance, function(x) {if(x){cols.TF[2]}else{cols.TF[1]}})
  #   barplot.list[[count]] = ggplot(fgseaResTidy, aes(reorder(pathway, psign), NES, fill = I(TFcolor))) +
  #     geom_col() + coord_flip() + theme_minimal() +
  #     labs(x="Pathway", y="Normalized Enrichment Score", title=cat)
  #   count = count + 1
  # }
  # print(wrap_plots(barplot.list, ncol = 2, nrow = 2))
  # ggsave(paste0(fig.dir,"/ps",i,"_GSEA_MSIGDB_",cancer,"_",j,"_DR.pdf"),
  #        width = 17, height = 10)
}



#find marker genes for each cluster
#nb "counts" is actually log counts
markers_default <- findMarkers(sce, 
                               groups = factor(sce$kmeans), # clusters to compare
                               test.type = "wilcox",   # t-test (default)
                               direction = "up", # test for either higher or lower expression (default)
                               lfc = 0, # null hypothesis log-fold-change = 0 (default)
                               pval.type = "all", # ranking of p-values based on any comparison (default)
                               assay.type = "counts"
)
#add cluster suffix to first three columns
markers_default = lapply(seq(length(markers_default)), function(x) {
  colnames(markers_default[[x]])[1:3] = paste(colnames(markers_default[[x]])[1:3],"_",
                                              names(markers_default)[x],sep="")
  markers_default[[x]] = as.data.frame(markers_default[[x]])
  markers_default[[x]]$gene = rownames(markers_default[[x]])
  return(markers_default[[x]][,c(1,2,3,ncol(markers_default[[x]]))])
})
markers_all = markers_default %>% reduce(left_join, by = "gene")
write.table(markers_all, file = paste0(fig.dir,"/cluster_marker_genes.tsv"), sep = "\t", 
            col.names = T, row.names = F, quote = F)

head(markers_all[order(markers_all$p.value_5, decreasing = F),],50)

#make trajectory figures with score defined states
df3 = df[sample(seq(nrow(df)),nrow(df)),]
p1 = ggplot(df3, aes(x = DR1, y = DR2, colour = sample)) + geom_point(size=pt.sz) + 
  guides(color = guide_legend(override.aes = list(size=2))) +
  theme(legend.key.size = unit(0.5, "cm"))
p2 = ggplot(df3, aes(x = DR1, y = DR2, colour = ICA_score_states)) + geom_point(size=pt.sz) + 
  guides(color = guide_legend(override.aes = list(size=2)))

cell.states = unique(df$ICA_score_states)
df.list = lapply(unique(df$sample), function(x) {
  data.frame(cell.state = cell.states,
             count = sapply(cell.states, function(y) {
               nrow(df[which(df$ICA_score_states==y & df$sample==x),])
             }),
             sample = x,
             grade = df$grade[which(df$sample==x)[1]])
})
df2 = do.call("rbind", df.list)
df2 = mutate(df2, sample = fct_reorder(sample, grade))
p3 = ggplot(df2, aes(fill=cell.state, y=count, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  theme(text=element_text(size=16)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))
plot_grid(p1, p2, p3, ncol = 3, align = "none")
ggsave(paste0(fig.dir,"/score_defined_cell_states.pdf"), height = 4, width = 18)

fwrite(df, file = paste0("df_file_traj_",cancer,".tsv.gz"))

###############################
#define cell state based on which branch cell belong to
###############################

#load df and sce from full dataset
load(paste0("post_traj_output_",cancer,"_all.RData"))

A = c(0, 0)
B = c(0, 1)
xs = c(A[1], B[1])
ys = c(A[2], B[2])
lm(ys ~ xs)

A = c(0, 0.005)
B = c(0.005,-0.005)
xs = c(A[1], B[1])
ys = c(A[2], B[2])
lm(ys ~ xs)

A = c(0, -0.007)
B = c(0.005, -0.003)
xs = c(A[1], B[1])
ys = c(A[2], B[2])
lm(ys ~ xs)

# df$cell.state.traj = sapply(seq(nrow(df)), function(x) {
#   if(df$DR2[x]> (0.005 - 2*df$DR1[x])){
#     "OPC-like"
#   }else if(df$DR1[x]< 0){
#     "astro-like"
#   }else if(df$DR2[x]< (-0.007 + 0.8*df$DR1[x])){
#     "neuro-like"
#   }else{
#     "intermediate"
#   }
# })
df$cell.state.traj = sapply(seq(nrow(df)), function(x) {
  if(df$DR2[x]> (0.0028 - 2*df$DR1[x])){
    "OPC-like"
  }else if(df$DR1[x]< (-0.0005)){
    "astro-like"
  }else if(df$DR2[x]< (-0.007 + 0.8*df$DR1[x])){
    "neuro-like"
  }else{
    "intermediate"
  }
})

fwrite(df, file = paste0("df_file_traj_",cancer,".tsv.gz"))

df3 = df[sample(seq(nrow(df)),nrow(df)),]
p1 = ggplot(df3, aes(x = DR1, y = DR2, colour = cell.state.traj)) + rasterise(geom_point(size=pt.sz),dpi=300) + 
  guides(color = guide_legend(override.aes = list(size=3))) + labs(colour = "Trajectory-defined\nmalignant states")+
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
p1 = add_curves(p1, inds)

cell.states = unique(df$cell.state.traj)
df.list = lapply(unique(df$sample), function(x) {
  data.frame(cell.state = cell.states,
             count = sapply(cell.states, function(y) {
               nrow(df[which(df$cell.state.traj==y & df$sample==x),])}),
             sample = x,
             grade = df$grade[which(df$sample==x)[1]])
})
df2 = do.call("rbind", df.list)
df2 = mutate(df2, sample = fct_reorder(sample, grade))
p2 = ggplot(df2, aes(fill=cell.state, y=count, x=sample)) + 
  geom_bar(position="fill", stat="identity") + ggnewscale::new_scale_fill() +
  geom_tile(aes(x = sample, y = 1.04, fill = factor(grade)), inherit.aes = FALSE, height = 0.06) +
  scale_fill_manual(name = "grade", values = c("2" = "skyblue","3" = "royalblue","4" = "darkblue"),
                    guide = guide_legend(order = 1)) +
  theme(text=element_text(size=16)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))
# p = (p1 | p2) + plot_layout(axes = "free")
# plot_grid(p1, p2, ncol = 2, align = "none")
ggsave(paste0(fig.dir,"/trajectory_defined_cell_states.pdf"), plot = p1, height = 4, width = 6)

df = df[,1:24]

#plot umap with cell types on it
df2 = df[sample(seq(nrow(df)),nrow(df)),]
ggplot(df2, aes(x = UMAP1, y = UMAP2, colour = cell.state.traj)) + geom_point(size=pt.sz) + 
  guides(color = guide_legend(override.aes = list(size=2)))
ggsave(paste0(fig.dir,"/umap_cell_states_traj_",cancer,".pdf"), width = 6, height = 4)


#use lamian approach to compare cell branch proportions by cell covariates, ini/rec, ini/rec paired, grade, treatment1

#compute proportion of each sample's cells in each cluster
branch.props = sapply(samp.md$Sample[which(!is.na(samp.md$Sample))], function(x) {
  sapply(unique(df$cell.state.traj), function(y) {
    return(sum(df$sample==x & df$cell.state.traj==y)/sum(df$sample==x))
  })
})
branch.props = branch.props[,grepl(cancer,colnames(branch.props))]
branch.props = branch.props[,apply(branch.props,2,function(x) {any(!is.na(x))})]


#ini rec
ini_rec.lam = sapply(colnames(branch.props), function(x) {unlist(str_split(x,"_" ))[3]})
id1 = which(ini_rec.lam=="ini")
id2 = which(ini_rec.lam=="rec1")
res.ir <- apply(branch.props, 1, function(i) t.test(i[id1], i[id2])$p.value)

#paired ini rec
branch.props = branch.props[,order(colnames(branch.props))]
pats = sapply(colnames(branch.props), function(x) {unlist(str_split(x,"_" ))[2]})
pats.paired = names(table(pats)[which(table(pats)==2)])
branch.props.paired = branch.props[,which(pats %in% pats.paired)]
res.ir.paired <- apply(branch.props.paired, 1, 
                       function(i) t.test(i[seq(1,ncol(branch.props.paired),2)], 
                                          i[seq(2,ncol(branch.props.paired),2)], paired=T)$p.value)

#grade
grade.lam = sapply(colnames(branch.props), function(x) {
  samp.md$Grade[which(samp.md$Sample==x)]
})
id1 = which(grade.lam==2)
id2 = which(grade.lam==3)
id3 = which(grade.lam==4)
res.gr23 <- apply(branch.props, 1, function(i) t.test(i[id1], i[id2])$p.value)
res.gr34 <- apply(branch.props, 1, function(i) t.test(i[id2], i[id3])$p.value)

#also do 2/3 v 4
id1 = which(grade.lam==2 | grade.lam==3)
id2 = which(grade.lam==4)
res.gr23_4 <- apply(branch.props, 1, function(i) t.test(i[id1], i[id2])$p.value)

branch.props.paired = branch.props.paired[order(rownames(branch.props.paired)),]
plot.list2 = list()
for(i in 1:nrow(branch.props.paired)){
  df.plot = data.frame(prop = branch.props.paired[i,], ini_rec = rep(c("initial","recurrence"),ncol(branch.props.paired)/2), 
                  patient = c(matrix(rep(pats.paired,2),ncol = length(pats.paired),byrow=T)))
  df.plot$ini_rec = factor(df.plot$ini_rec, levels=c("initial","recurrence"))
  df.plot$patient = factor(df.plot$patient)
  #Wilcoxon test
  pval = wilcox.test(df.plot$prop[seq(1,nrow(df.plot),2)], df.plot$prop[seq(2,nrow(df.plot),2)], paired = TRUE)$p.value
  
  plot.list2[[i]] = ggplot(df.plot, aes(x = ini_rec, y = prop)) + geom_boxplot(outlier.shape=NA,fill = gg_color_hue(nrow(branch.props.paired))[i]) +
    geom_line(aes(group = patient)) +
    geom_point(size = 2) + labs(y="Proportion",x="") + theme(text = element_text(size = 19)) +
    ggtitle(paste0(rownames(branch.props.paired)[i]," (p=", signif(pval, digits=2),")"))
}
p1 = wrap_plots(plot.list2)
print(p1)
ggsave(paste0(fig.dir,"/branch_proportions_ini_rec_paired.pdf"), width = 8, height = 8)


comparisons <- list(c("2","3"), c("3","4"), c("2","4"))

df.plot = do.call("rbind", lapply(seq(nrow(branch.props)), function(i) {
  data.frame(prop = branch.props[i,], grade = grade.lam, cell_state = rownames(branch.props)[i],
             grade_comb = sapply(grade.lam, function(x) {if(x==2 || x==3){"2/3"}else if(x==4){"4"}}))
}))
df.plot$grade = factor(df.plot$grade)
# ggplot(df.plot, aes(x = grade, y = prop, fill = cell_state)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.8) +
#   geom_jitter(width = 0.15, size = 1.5) +
#   # Wilcoxon test
#   stat_compare_means(comparisons = comparisons,method = "wilcox.test",label = "p.signif", group.by = "cell_state") +
#   facet_wrap(~cell_state, ncol=4) +   # makes the 4 panels like your figure
#   theme_classic() + labs(x = "Grade",y = "Proportion in branch")
# ggsave(paste0(fig.dir,"/branch_proportions_grade_not_comb.pdf"), width = 8, height = 3)

library(rstatix)

stats <- df.plot %>%
  group_by(cell_state) %>%
  wilcox_test(prop ~ grade, comparisons = comparisons) %>%
  add_significance(p.col = "p") %>%     # this creates p.signif
  add_xy_position(x = "grade")

ggplot(df.plot, aes(grade, prop, fill = cell_state)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width = 0.15) +
  stat_pvalue_manual(stats, label = "p.signif") +
  facet_wrap(~cell_state, ncol = 4) +
  theme_classic() + labs(y="Proportion")
ggsave(paste0(fig.dir,"/branch_proportions_grade_not_comb.pdf"), width = 8, height = 3)




########################################
##scenicplus analysis
########################################


load(paste0("post_traj_output_",cancer,"_",j,rem,".RData"))
df = fread(file = paste0("df_file_traj_",cancer,".tsv.gz"), header = T, data.table=F)
rownames(df) = df$cell

out.dir = paste0("/Users/alyne/Library/CloudStorage/Dropbox/Curie/Curie_StCloud/anne-marie/10x_RNA_ATAC_oligo_analysis",
                 "/SCENICplus_soroban/",cancer,"_full_run/figures_scenicplus_newdb_",cancer)

use.extended = T
#read in region AUC scores
TF.AUCr.d = fread(paste0(out.dir,"/direct_region_based_AUC_",cancer,".tsv.gz"), header=T, data.table = F)
if(use.extended){
  TF.AUCr.e = fread(paste0(out.dir,"/extended_region_based_AUC_",cancer,".tsv.gz"), header=T, data.table = F)
}else{
  TF.AUCr.e = data.frame(Cell = TF.AUCr.d$Cell)
}
TF.AUCr = cbind(TF.AUCr.d,TF.AUCr.e[,-1])
rownames(TF.AUCr) = TF.AUCr$Cell
rownames(TF.AUCr) = sapply(rownames(TF.AUCr), function(x) {unlist(strsplit(x,"___",fixed=T))[1]})
TF.AUCr$Cell = NULL


#some summary figures of AUC scores and numbers of genes, types of regulons etc
type.reg = sapply(colnames(TF.AUCr), function(x) {unlist(strsplit(x,"_"))[3]})
df.types = data.frame(type = names(table(type.reg)), count = table(type.reg))
ggplot(df.types, aes(x = type, y = count.Freq)) + geom_bar(stat="identity")
ggsave(paste0(fig.dir,"/barplot_regulon_types_",cancer,".pdf"))

#remove regulons with negative TF region
to.rem = type.reg %in% c("+/-","-/-")

#find cells in both dfs
joint.cells = intersect(rownames(df),rownames(TF.AUCr))
df.ss.reg = df[sapply(joint.cells, function(x) {which(rownames(df)==x)}),]
TF.AUCr.malig = TF.AUCr[sapply(joint.cells, function(x) {which(rownames(TF.AUCr)==x)}),]
TF.AUCr.malig = TF.AUCr.malig[,!to.rem]

#for both pseudotimes
for(i in 1:2){
  #compute corelation
  corr_ps = lapply(seq(ncol(TF.AUCr.malig)), function(x) {
    res = cor.test(df.ss.reg[,which(colnames(df.ss.reg)==paste0("slingPseudotime_",i))],TF.AUCr.malig[,x], 
                   use = "complete.obs", method = "pearson")
    return(c(res$estimate, res$p.value))
  })
  
  df.corr = do.call("rbind",corr_ps)
  df.corr = as.data.frame(df.corr)
  rownames(df.corr) = colnames(TF.AUCr.malig)
  colnames(df.corr) = c("corr_ps","pval_1")
  
  Nplot = 8
  for(direc in c("up","down")){
    plot.list = list()
    plot.list2 = list()
    #plot top 5, both expression on trajectory and dot plot of expr v pseudotime
    for(b in 1:Nplot){
      if(direc=="up"){
        gene = rownames(df.corr[order(df.corr$corr_ps,decreasing=T),])[b]
      }else{
        gene = rownames(df.corr[order(df.corr$corr_ps),])[b]
      }
      df3 = data.frame(pseudotime = df.ss.reg[,which(colnames(df.ss.reg)==paste0("slingPseudotime_",i))], 
                       activity = TF.AUCr.malig[,which(colnames(TF.AUCr.malig)==gene)],
                       cluster = df.ss.reg$kmeans_DR)
      plot.list2[[b]] = ggplot(df3, aes(x = pseudotime, y = activity, color = cluster)) +
        geom_point(size=pt.sz) + labs(title=gene)
      
      df3 = cbind(df3, DR1=df.ss.reg$DR1, DR2=df.ss.reg$DR2)
      plot.list[[b]] = ggplot(df3, aes(x = DR1, y = DR2, color = activity)) + geom_point(size=pt.sz) #+
      # geom_path(data = df.crv, aes(x=crv1_1,y=crv1_2),col="black") +
      # geom_path(data = df.crv, aes(x=crv2_1,y=crv2_2),col="black")
    }
    p1 = wrap_plots(c(plot.list,plot.list2), ncol = Nplot, byrow = T)
    print(p1)
    ggsave(paste0(fig.dir,"/ps",i,"_TF_reg_",direc,"_",cancer,"_",j,"_",dr,".pdf"),
           width = 25, height = 5)
  }
}


#read in gene AUC scenicplus output
TF.AUCg.d = fread(paste0(out.dir,"/direct_gene_based_AUC_",cancer,".tsv.gz"), header=T, data.table = F)
if(use.extended){
  TF.AUCg.e = fread(paste0(out.dir,"/extended_gene_based_AUC_",cancer,".tsv.gz"), header=T, data.table = F)
}else{
  TF.AUCg.e = data.frame(Cell = TF.AUCg.d$Cell)
}
TF.AUCg = cbind(TF.AUCg.d,TF.AUCg.e[,-1])
rownames(TF.AUCg) = TF.AUCg$Cell
rownames(TF.AUCg) = sapply(rownames(TF.AUCg), function(x) {unlist(strsplit(x,"___",fixed=T))[1]})
TF.AUCg$Cell = NULL
#find cells in both dfs
joint.cells = intersect(rownames(df),rownames(TF.AUCg))
df.ss = df[sapply(joint.cells, function(x) {which(rownames(df)==x)}),]
TF.AUCg.malig = TF.AUCg[sapply(joint.cells, function(x) {which(rownames(TF.AUCg)==x)}),]
TF.AUCg.malig = TF.AUCg.malig[,!to.rem]


Ngene.reg = sapply(colnames(TF.AUCg), function(x) {
  str = unlist(strsplit(x,"_"))[4]
  as.numeric(substr(str,2,nchar(str)-2))
})
Nreg.reg = sapply(colnames(TF.AUCr), function(x) {
  str = unlist(strsplit(x,"_"))[4]
  as.numeric(substr(str,2,nchar(str)-2))
})
df.ngene = data.frame(Ngene = Ngene.reg, Nregion = Nreg.reg)
ggplot(df.ngene, aes(x = Ngene)) + geom_histogram()
ggsave(paste0(fig.dir,"/barplot_regulon_Ngene_",cancer,".pdf"))
ggplot(df.ngene, aes(x = Nregion)) + geom_histogram()
ggsave(paste0(fig.dir,"/barplot_regulon_Nregion_",cancer,".pdf"))
ggplot(df.ngene, aes(x = Nregion, y = Ngene)) + geom_point()
ggsave(paste0(fig.dir,"/scatter_regulon_Nregion_v_Ngene_",cancer,".pdf"))

#AUC scores per gene/region set size
df.AUC = do.call("rbind",lapply(seq(ncol(TF.AUCr)), 
                                function(x) {data.frame(AUCr = TF.AUCr[,x], Nregion = Nreg.reg[x])}))
df.AUC$Nreg_group <- cut(df.AUC$Nregion,
                         breaks = c(1,500,1000,1500,2000,2500,11000),
                         include.lowest = TRUE)
ggplot(df.AUC, aes(x = Nreg_group, y = AUCr)) +
  geom_boxplot() + xlab("Nreg bins") + ylab("AUCr") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(fig.dir,"/boxplot_AUC_region_scores_",cancer,".pdf"))

df.AUC = do.call("rbind",lapply(seq(ncol(TF.AUCg)), 
                                function(x) {data.frame(AUCg = TF.AUCg[,x], Ngene = Ngene.reg[x])}))
df.AUC$Ngene_group <- cut(df.AUC$Ngene,
                          breaks = c(1,200,400,600,800,1000,3000),
                          include.lowest = TRUE)
ggplot(df.AUC, aes(x = Ngene_group, y = AUCg)) +
  geom_boxplot() + xlab("Ngene bins") + ylab("AUCg") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(fig.dir,"/boxplot_AUC_genes_scores_",cancer,".pdf"))


#correlation with pseudotime
df.corr = list()
for(i in 1:2){
  #compute corelation
  corr_ps = lapply(seq(ncol(TF.AUCg.malig)), function(x) {
    res = cor.test(df.ss[,which(colnames(df.ss)==paste0("slingPseudotime_",i))],TF.AUCg.malig[,x], 
                   use = "complete.obs", method = "pearson")
    return(c(res$estimate, res$p.value))
  })
  
  df.corr[[i]] = do.call("rbind",corr_ps)
  df.corr[[i]] = as.data.frame(df.corr[[i]])
  rownames(df.corr[[i]]) = colnames(TF.AUCg.malig)
  colnames(df.corr[[i]]) = c("corr_ps","pval_1")
  
  Nplot = 6
  for(direc in c("up","down")){
    plot.list = list()
    plot.list2 = list()
    #plot top 5, both expression on trajectory and dot plot of expr v pseudotime
    for(b in 1:Nplot){
      if(direc=="up"){
        gene = rownames(df.corr[[i]][order(df.corr[[i]]$corr_ps,decreasing=T),])[b]
      }else{
        gene = rownames(df.corr[[i]][order(df.corr[[i]]$corr_ps),])[b]
      }
      df3 = data.frame(pseudotime = df.ss[,which(colnames(df.ss)==paste0("slingPseudotime_",i))], 
                       activity = TF.AUCg.malig[,which(colnames(TF.AUCg.malig)==gene)],
                       cluster = df.ss$kmeans_DR)
      plot.list2[[b]] = ggplot(df3, aes(x = pseudotime, y = activity, color = cluster)) +
        geom_point(size=pt.sz) + labs(title=gene)
      
      df3 = cbind(df3, DR1=df.ss$DR1, DR2=df.ss$DR2)
      plot.list[[b]] = ggplot(df3, aes(x = DR1, y = DR2, color = activity)) + geom_point(size=pt.sz) #+
      # geom_path(data = df.crv, aes(x=crv1_1,y=crv1_2),col="black") +
      # geom_path(data = df.crv, aes(x=crv2_1,y=crv2_2),col="black")
    }
    p1 = wrap_plots(c(plot.list,plot.list2), ncol = Nplot, byrow = T)
    print(p1)
    ggsave(paste0(fig.dir,"/ps",i,"_TF_act_",direc,"_",cancer,"_",j,"_",dr,".pdf"),
           width = 20, height = 5)
  }
}


states = c("neuro-like","astro-like","OPC-like","intermediate")

state.meds = list()
oneway.list = list()
df.anova = list()
for(i in 1:ncol(TF.AUCg.malig)){
  df.anova[[i]] = data.frame(cell_state = df.ss$cell.state.traj, auc = TF.AUCg.malig[,i])
  state.meds[[i]] = sapply(states, function(x) {
    median(df.anova[[i]]$auc[which(df.anova[[i]]$cell_state==x)])
  })
  oneway.list[[i]] = oneway.test(auc ~ cell_state, data = df.anova[[i]])
}
names(state.meds) = colnames(TF.AUCg.malig)

maxs = sapply(state.meds, function(x) {which.max(x)})
#all signifcant 
#filter to take only those with neuro median highest
#then order by biggest difference with second largest
for(cell_state in c(1,2,3,4)){
  state.meds.filt = state.meds[which(maxs==cell_state)]
  oneway.list.filt = oneway.list[which(maxs==cell_state)]
  df.anova.filt = df.anova[which(maxs==cell_state)]
  diffs = sapply(state.meds.filt, function(x) {min(x[cell_state] - x[c(1,2,3,4)[-cell_state]])})
  inds = sort.int(diffs,decreasing = T,index.return = T)
  
  for(i in 1:length(inds$ix[1:20])){
    df.to.plot = df.anova.filt[[inds$ix[i]]]
    df.to.plot = do.call("rbind",lapply(states, function(x) {
      df = df.to.plot[df.to.plot$cell_state==x,]
      if(nrow(df)>3000){
        df = df[sample(nrow(df),3000),]
      }
      return(df)
    }))
    ggplot(df.to.plot, aes(x = cell_state, y = auc)) + geom_boxplot() + geom_jitter(size=0.01) + 
      labs(title=names(state.meds.filt)[inds$ix[i]])
    ggsave(paste0(fig.dir,"/",states[cell_state],"_anova_TF_gene_",i,".pdf"))
  }
}



#try a different approach using a linear model so can include patient and sample effects
library(lme4)
library(lmerTest)


state.meds = list()
lmer.fit = list()
df.anova = list()
for(i in 1:ncol(TF.AUCg.malig)){
  df.anova[[i]] = data.frame(cell_state = df.ss$cell.state.traj, auc = TF.AUCg.malig[,i], 
                             sample = df.ss$sample, patient = sapply(df.ss$sample, function(x) {unlist(strsplit(x,"_"))[2]}),
                             timepoint = df.ss$ini_rec)
  state.meds[[i]] = sapply(states, function(x) {
    median(df.anova[[i]]$auc[which(df.anova[[i]]$cell_state==x)])
  })
  lmer.fit[[i]] <- lmer(auc ~ cell_state + timepoint + (1|sample), data = df.anova[[i]])
}
names(state.meds) = colnames(TF.AUCg.malig)
names(lmer.fit) = colnames(TF.AUCg.malig)

#for each state, extract genes with biggest up change in each state
fixed.effects = lapply(lmer.fit, function(x) {
  c(fixef(x)[2:4],"cell_stateastro-like"=0)
})
names(fixed.effects) = names(lmer.fit)

maxs = sapply(fixed.effects, function(x) {which.max(x)})
for(cell_state in c(1,2,3,4)){
  fixed.effects.filt = fixed.effects[which(maxs==cell_state)]
  df.anova.filt = df.anova[which(maxs==cell_state)]
  diffs = sapply(fixed.effects.filt, function(x) {min(x[cell_state] - x[c(1,2,3,4)[-cell_state]])})
  inds = sort.int(diffs,decreasing = T,index.return = T)
  
  for(i in 1:length(inds$ix[1:20])){
    df.to.plot = df.anova.filt[[inds$ix[i]]]
    df.to.plot = do.call("rbind",lapply(states, function(x) {
      df = df.to.plot[df.to.plot$cell_state==x,]
      if(nrow(df)>3000){
        df = df[sample(nrow(df),3000),]
      }
      return(df)
    }))
    ggplot(df.to.plot, aes(x = cell_state, y = auc)) + geom_boxplot() + geom_jitter(size=0.01) + 
      labs(title=names(fixed.effects.filt)[inds$ix[i]])
    ggsave(paste0(fig.dir,"/",substr(names(fixed.effects[[1]])[cell_state],11,nchar(names(fixed.effects[[1]])[cell_state])),
                  "_lmer_TF_gene_",i,".pdf"))
  }
}


#also look at TFs with largest timepoint impact
fixed.effects.time = sapply(lmer.fit, function(x) {
  fixef(x)[5]
})
names(fixed.effects.time) = names(lmer.fit)
fixed.effects.time = sort(fixed.effects.time, decreasing=T)
for(i in 1:20){
  df.anova = data.frame(cell_state = df.ss$cell.state.traj, sample = df.ss$sample, timepoint = df.ss$ini_rec,
                        auc = TF.AUCg.malig[,which(colnames(TF.AUCg.malig)==names(fixed.effects.time)[i])])
  df.to.plot = do.call("rbind",lapply(states, function(x) {
    df = df.anova[df.anova$cell_state==x,]
    if(nrow(df)>3000){
      df = df[sample(nrow(df),3000),]
    }
    return(df)
  }))
  ggplot(df.to.plot, aes(x = cell_state, y = auc, fill = timepoint)) + geom_boxplot(position = position_dodge(width = 0.75)) +
    geom_jitter(size=0.01,position = position_jitterdodge(jitter.width = 0.15,dodge.width = 0.75)) +
    labs(title=names(fixed.effects.time)[inds$ix[i]])
  ggsave(paste0(fig.dir,"/timepoint_lmer_TF_up_gene_",i,".pdf"), width=6)
}

fixed.effects.time = sort(fixed.effects.time, decreasing=F)
for(i in 1:20){
  df.anova = data.frame(cell_state = df.ss$cell.state.traj, sample = df.ss$sample, timepoint = df.ss$ini_rec,
                        auc = TF.AUCg.malig[,which(colnames(TF.AUCg.malig)==names(fixed.effects.time)[i])])
  df.to.plot = do.call("rbind",lapply(states, function(x) {
    df = df.anova[df.anova$cell_state==x,]
    if(nrow(df)>3000){
      df = df[sample(nrow(df),3000),]
    }
    return(df)
  }))
  ggplot(df.to.plot, aes(x = cell_state, y = auc, fill = timepoint)) + geom_boxplot(position = position_dodge(width = 0.75)) +
    geom_jitter(size=0.01,position = position_jitterdodge(jitter.width = 0.15,dodge.width = 0.75)) +
    labs(title=names(fixed.effects.time)[inds$ix[i]])
  ggsave(paste0(fig.dir,"/timepoint_lmer_TF_down_gene_",i,".pdf"), width=6)
}


#plot transcription factors as heatmap along pseudotime
library(colorRamp2)
library(ComplexHeatmap)


################
## tradeseq analysis
################
library(tradeSeq)
library(tidyr)

i=2
pseudotime_mat <- matrix(df.ss[,which(colnames(df.ss)==paste0("slingPseudotime_",i))], ncol = 1)
indsNA = which(!is.na(pseudotime_mat))
cellWeights <- matrix(1, nrow = nrow(pseudotime_mat), ncol = 1)

gam_fit <- fitGAM(counts = t(TF.AUCg.malig[indsNA,]), pseudotime = pseudotime_mat[indsNA,], 
                  cellWeights = cellWeights[indsNA,], nknots = 6)
smooth_vals <- predictSmooth(gam_fit, colnames(TF.AUCg.malig))
assoc <- associationTest(gam_fit)

#find significant regulons
sig_regs <- rownames(assoc)[assoc$pvalue < 0.005 & assoc$meanLogFC > 0]
smooth_sig <- smooth_vals[smooth_vals$gene %in% sig_regs,]
#convert to matrix
heat_df <- smooth_sig %>%
  select(gene, time, yhat) %>%
  pivot_wider(names_from = time, values_from = yhat)
heat_mat <- as.matrix(heat_df[,-1])
rownames(heat_mat) <- heat_df$gene

#remove those where maximum value is below 0.015
peak_idx <- apply(heat_mat, 1, which.max)
to.keep = which(heat_mat[cbind(seq(nrow(heat_mat)),peak_idx)]>0.015)
peak_idx = peak_idx[to.keep]
heat_mat = heat_mat[to.keep,]

col_fun <- colorRamp2(c(0, 0.05, 0.2), c("lightgrey", "orange","blue"))
ht = Heatmap(
  heat_mat,
  name = "AUC",
  cluster_rows = T,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  col = col_fun,
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
pdf(paste0(fig.dir,"/ps",i,"_heatmap_TFactivity_tradeseq_",cancer,".pdf"), width = 15, height = 12, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


col_fun <- colorRamp2(c(0, 0.5, 1), c("lightgrey", "orange","blue"))
#same but scaled
heat_mat2 = t(apply(heat_mat, 1, function(x)(x-min(x))/(max(x)-min(x))))
ht = Heatmap(
  heat_mat2,
  name = "AUC",
  cluster_rows = T,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  col = col_fun,
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
pdf(paste0(fig.dir,"/ps",i,"_heatmap_TFactivity_tradeseq_scaled_",cancer,".pdf"), width = 15, height = 12, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


#put into three clusters, peak in middle, up or down
#first identify and extract those that peak in middle
pt <- as.numeric(colnames(heat_mat))
peak_pt <- pt[peak_idx]
mid_regs <- names(peak_idx)[
  peak_pt > quantile(pt, 0.35) &
    peak_pt < quantile(pt, 0.65)
]
heat_mat.peak = heat_mat[which(rownames(heat_mat) %in% mid_regs),]
#split the rest into up or down
heat_mat.rest = heat_mat[which(!(rownames(heat_mat) %in% mid_regs)),]
heat_mat.up = heat_mat.rest[heat_mat.rest[,1] < heat_mat.rest[,ncol(heat_mat.rest)],]
heat_mat.down = heat_mat.rest[heat_mat.rest[,1] >= heat_mat.rest[,ncol(heat_mat.rest)],]
heat_mat.split = rbind(heat_mat.peak,heat_mat.up,heat_mat.down)

row_anno_df = data.frame(regulon_type = c(rep("peaked",nrow(heat_mat.peak)),
                                          rep("up",nrow(heat_mat.up)),
                                          rep("down",nrow(heat_mat.down))))
rownames(row_anno_df) = rownames(heat_mat.split)
left_anno <- rowAnnotation(df = row_anno_df)


col_fun <- colorRamp2(c(0, 0.05, 0.2), c("lightgrey", "orange","blue"))
ht = Heatmap(
  heat_mat.split,
  name = "AUC",
  # cluster_rows = F,
  row_split = row_anno_df$regulon_type,
  # cluster_row_slices = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  left_annotation = left_anno,
  col = col_fun,
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
pdf(paste0(fig.dir,"/ps",i,"_heatmap_TFactivity_tradeseq_split_",cancer,".pdf"), width = 15, height = 12, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


#smaller plot with just 5 for each category
Npertype = 8
Nuse=20
Ncol = ncol(heat_mat)
#split the rest into up or down
heat_mat.scale = t(apply(heat_mat, 1, function(x)(x-min(x))/(max(x)-min(x))))
rownames(heat_mat.scale) = rownames(heat_mat)
heat_mat.peak = heat_mat.scale[which(rownames(heat_mat.scale) %in% mid_regs),]
heat_mat.rest = heat_mat.scale[which(!(rownames(heat_mat.scale) %in% mid_regs)),]
heat_mat.up = heat_mat.rest[rowMeans(heat_mat.rest[,1:Nuse]) < rowMeans(heat_mat.rest[,(ncol(heat_mat.rest)-Nuse+1):ncol(heat_mat.rest)]),]
heat_mat.down = heat_mat.rest[rowMeans(heat_mat.rest[,1:Nuse]) >= rowMeans(heat_mat.rest[,(ncol(heat_mat.rest)-Nuse+1):ncol(heat_mat.rest)]),]
heat_mat.split = rbind(heat_mat.peak[order(rowMeans(cbind(rowMeans(heat_mat.peak[,40:60])-rowMeans(heat_mat.peak[,1:Nuse]),rowMeans(heat_mat.peak[,40:60])-rowMeans(heat_mat.peak[,(Ncol-Nuse+1):Ncol]))),
                                           decreasing=T)[1:Npertype],],
                       heat_mat.up[order(rowMeans(heat_mat.peak[,(Ncol-Nuse+1):Ncol])-rowMeans(heat_mat.peak[,1:Nuse]),decreasing=T)[1:Npertype],],
                       heat_mat.down[order(rowMeans(heat_mat.peak[,1:Nuse])-rowMeans(heat_mat.peak[,(Ncol-Nuse+1):Ncol]),decreasing=T)[1:Npertype],])
row_anno_df = data.frame(regulon_type = c(rep("peaked",Npertype), rep("up",Npertype), rep("down",Npertype)))
rownames(row_anno_df) = rownames(heat_mat.split)
left_anno <- rowAnnotation(df = row_anno_df)
col_fun <- colorRamp2(c(0, 0.5, 1), c("lightgrey", "orange","blue"))
heat_mat.split2 = t(apply(heat_mat.split, 1, function(x)(x-min(x))/(max(x)-min(x))))
rownames(heat_mat.split2) = sapply(rownames(heat_mat.split2), function(x) {unlist(  strsplit(x,"_")  )[1]})
ht = Heatmap(
  heat_mat.split2,
  name = "AUC",
  # cluster_rows = F,
  row_split = row_anno_df$regulon_type,
  # cluster_row_slices = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  left_annotation = left_anno,
  col = col_fun,
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
pdf(paste0(fig.dir,"/ps",i,"_heatmap_TFactivity_tradeseq_split_small_scaled_",cancer,".pdf"), width = 10, height = 5, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


#fit gam so can use gaussian rather than nb

library(mgcv)

gam_fit.list = list()
for(kk in 1:ncol(TF.AUCg.malig)){
  df = data.frame(auc = TF.AUCg.malig[indsNA,kk], pseudotime = pseudotime_mat[indsNA,],
                  sample = df.ss$sample[indsNA])
  df$sample = factor(df$sample)
  df$auc = df$auc - min(df$auc)
  df$auc = df$auc * (1/(max(df$auc)))
  gam_fit.list[[kk]] = gam(auc ~ s(pseudotime, k = 3) + s(sample, bs = "re"),
    data = df, family = gaussian()
  )
}


res = data.frame()
for(kk in seq_along(gam_fit.list)) {
  s.tab <- summary(gam_fit.list[[kk]])$s.table
  res <- rbind(res,
    data.frame(TF = colnames(TF.AUCg.malig)[kk], edf = s.tab["s(pseudotime)", "edf"],
      F = s.tab["s(pseudotime)", "F"], pval = s.tab["s(pseudotime)", "p-value"],
      r2 = summary(gam_fit.list[[kk]])$r.sq))
}
res$padj <- p.adjust(res$pval, method = "BH")


###need to sort by pseudotime
ps.order =  sort.int(pseudotime_mat[indsNA,], index.return = T)
res$delta = sapply(gam_fit.list, function(x) {
  mean(x$fitted.values[ps.order$ix[(length(ps.order$ix)-9):length(ps.order$ix)]]) - 
    mean(x$fitted.values[ps.order$ix[1:10]])
})
sig <- subset(res,padj < 0.005)
sig = sig[order(sig$delta,decreasing = T),]

Nplot = 20
# inds = match(sig$TF[c(1:Nplot,(length(sig$TF)-Nplot+1):length(sig$TF))],colnames(TF.AUCg.malig))
inds = match(sig$TF[c(1:Nplot)],colnames(TF.AUCg.malig))
heat_mat.gam = do.call("rbind",lapply(gam_fit.list[inds], function(x) {x$fitted.values[ps.order$ix]}))
rownames(heat_mat.gam) = sig$TF[c(1:Nplot)]
ps.ordered = pseudotime_mat[indsNA,][ps.order$ix]

#get cell indices to sample evenly along pseudotime
Nsamp = 1000
times = seq(0,max(pseudotime_mat[indsNA,]),max(pseudotime_mat[indsNA,])/Nsamp)
inds.ord = sapply(times, function(x) {which.min(abs(ps.ordered-x))})
column_ha = HeatmapAnnotation(cell_state = df.ss$cell.state.traj[indsNA][ps.order$ix][inds.ord])

col_fun <- colorRamp2(c(0, 0.3, 0.8), c("lightgrey", "orange","blue"))
ht = Heatmap(
  heat_mat.gam[,inds.ord],
  name = "AUC",
  cluster_rows = T,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  bottom_annotation = column_ha,
  col = col_fun,
  rect_gp = gpar(col = NA),
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y,
      w + unit(0.3, "pt"),
      h + unit(0.3, "pt"),
      gp = gpar(fill = fill, col = NA))
  },use_raster = FALSE
)
pdf(paste0(fig.dir,"/ps",i,"_heatmap_TFactivity_gam_",cancer,".pdf"), width = 15, height = 6, useDingbats = FALSE)
ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()




# #try identifying using wilcox test
# groups <- df$cell.state.traj
# unique_groups <- unique(groups)
# group_indices <- lapply(unique_groups, function(g) which(groups == g))
# names(group_indices) <- unique_groups
# 
# # Main loop
# mat = as.matrix(TF.AUCg.malig)
# regs.spec <- lapply(unique_groups, function(g) {
#   idx_in  <- group_indices[[g]]
#   idx_out <- which(groups != g)
#   
#   pvals <- vapply(seq_len(ncol(TF.AUCg.malig)), function(i) {
#     y <- mat[,i]
#     out <- try(
#       wilcox.test(y[idx_in], y[idx_out], alternative = "less", exact = FALSE, correct = FALSE), silent = TRUE
#     )
#     
#     if (inherits(out, "try-error")) return(NA_real_)
#     out$p.value
#   }, numeric(1))
#   
#   colnames(TF.AUCg.malig)[order(pvals)[1:Npertype]]
# })
# 
# smooth_sig <- smooth_vals[smooth_vals$gene %in% unlist(regs.spec),]
# #convert to matrix
# heat_df <- smooth_sig %>%
#   select(gene, time, yhat) %>%
#   pivot_wider(names_from = time, values_from = yhat)
# heat_mat <- as.matrix(heat_df[,-1])
# rownames(heat_mat) <- heat_df$gene
# heat_mat.split = do.call("rbind",lapply(seq(length(regs.spec)), function(x) {
#   rel.rows = heat_mat[regs.spec[[x]],]
#   rownames(rel.rows) = paste0(rownames(rel.rows),"-",x)
#   return(rel.rows)
#   }))
# # heat_mat.split = do.call("rbind",lapply(regs.spec, function(x) {heat_mat[x %in% rownames(heat_mat),]}))
# 
# row_anno_df = data.frame(regulon_type = unlist(lapply(unique_groups,function(x) {rep(x,Npertype)})))
# rownames(row_anno_df) = rownames(heat_mat.split)
# left_anno <- rowAnnotation(df = row_anno_df)
# 
# 
# col_fun <- colorRamp2(c(0, 0.05, 0.2), c("lightgrey", "orange","blue"))
# ht = Heatmap(
#   heat_mat.split,
#   name = "AUC",
#   # cluster_rows = F,
#   row_split = row_anno_df$regulon_type,
#   # cluster_row_slices = TRUE,
#   cluster_columns = FALSE,
#   show_row_names = TRUE,
#   show_column_names = FALSE,
#   row_names_gp = gpar(fontsize = 8),
#   left_annotation = left_anno,
#   col = col_fun,
#   rect_gp = gpar(col = NA),
#   layer_fun = function(j, i, x, y, w, h, fill) {
#     grid.rect(
#       x, y,
#       w + unit(0.3, "pt"),
#       h + unit(0.3, "pt"),
#       gp = gpar(fill = fill, col = NA)
#     )
#   }
# )
# pdf(paste0(fig.dir,"/ps",i,"_heatmap_TFactivity_tradeseq_split_small_",cancer,".pdf"), width = 10, height = 5, useDingbats = FALSE)
# draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
# dev.off()
# 
# #same but scaled
# col_fun <- colorRamp2(c(0, 0.5, 1), c("lightgrey", "orange","blue"))
# heat_mat.split2 = t(apply(heat_mat.split, 1, function(x)(x-min(x))/(max(x)-min(x))))
# rownames(heat_mat.split2) = sapply(rownames(heat_mat.split2), function(x) {unlist(  strsplit(x,"_")  )[1]})
# ht = Heatmap(
#   heat_mat.split2,
#   name = "AUC",
#   # cluster_rows = F,
#   row_split = row_anno_df$regulon_type,
#   # cluster_row_slices = TRUE,
#   cluster_columns = FALSE,
#   show_row_names = TRUE,
#   show_column_names = FALSE,
#   row_names_gp = gpar(fontsize = 8),
#   left_annotation = left_anno,
#   col = col_fun,
#   rect_gp = gpar(col = NA),
#   layer_fun = function(j, i, x, y, w, h, fill) {
#     grid.rect(
#       x, y,
#       w + unit(0.3, "pt"),
#       h + unit(0.3, "pt"),
#       gp = gpar(fill = fill, col = NA)
#     )
#   }
# )
# pdf(paste0(fig.dir,"/ps",i,"_heatmap_TFactivity_tradeseq_split_small_scaled_",cancer,".pdf"), width = 10, height = 5, useDingbats = FALSE)
# draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
# dev.off()









#fit gam for regions
gam_fit.reg <- fitGAM(counts = t(TF.AUCr.malig[indsNA,]), pseudotime = pseudotime_mat[indsNA,], 
                  cellWeights = cellWeights[indsNA,], nknots = 6)
smooth_vals.reg <- predictSmooth(gam_fit.reg, colnames(TF.AUCr.malig))
assoc <- associationTest(gam_fit.reg)

#normalise
scale01 <- function(x) (x - min(x)) / (max(x) - min(x))

auc_n  <- scale01(smooth_vals.reg)
expr_n <- scale01(smooth_vals)

#compute lag using cross correlation
ccf_res <- ccf(auc_n, expr_n, plot = FALSE)

lag <- ccf_res$lag[which.max(ccf_res$acf)]

# lag > 0 → ATAC leads RNA
# regions open before gene expression
# lag < 0 → RNA leads ATAC (rare biologically)
# lag ≈ 0 → simultaneous

#convert lag to same units as pseudotime
pt_lag <- lag * (diff(range(pt_grid)) / length(pt_grid))

#loop through transcription factors
lags <- sapply(pairs, function(pair) {
  auc_pred  <- predict(pair$gam_auc,  newdata = data.frame(pt = pt_grid))
  expr_pred <- predict(pair$gam_expr, newdata = data.frame(pt = pt_grid))
  
  auc_n  <- scale01(auc_pred)
  expr_n <- scale01(expr_pred)
  
  ccf_res <- ccf(auc_n, expr_n, plot = FALSE)
  ccf_res$lag[which.max(ccf_res$acf)]
})









#for each pseudotime
for(i in 1:2){
  #subsample along trajectory and reorder
  mx.time = max(df.ss[[paste0("slingPseudotime_",i)]],na.rm=T)
  Nss = 2000
  ts = seq(0,mx.time,length.out=Nss)
  cells = unique(sapply(ts, function(x) {df.ss$cell[which.min(abs(df.ss[[paste0("slingPseudotime_",i)]]-x))]}))
  df.toplot = TF.AUCg[sapply(cells, function(x) {which(rownames(TF.AUCg)==x)}),]
  
  col_fun <- colorRamp2(c(0, 0.1, 0.25), c("lightgrey", "orange","blue"))
  ht <- Heatmap(
    t(df.toplot),
    name = "TFactivity",
    col = col_fun,
    cluster_rows = T,
    cluster_columns = F,
    row_gap = unit(0, "mm"),
    column_gap = unit(0, "mm"),
    show_row_names = T,
    row_names_gp = gpar(fontsize = 3),
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
  pdf(paste0(fig.dir,"/ps",i,"_heatmap_TFactivity_",cancer,".pdf"), width = 15, height = 20, useDingbats = FALSE)
  draw(ht)
  dev.off()
}



#read in rss scores
rss = fread(paste0(out.dir,"/figures/rss_",cancer,".csv"), header=T, data.table = F)
rss = rss[,c(1,c(2:ncol(rss))[!to.rem])]

#get top positively and negatively correlated, plus top of each state by RSS, top by activity 
#and known TFs
top.rss = unlist(lapply(c("OPC-like","astro-like","intermediate","neuro-like"), function(x) {
  y = as.numeric(rss[which(rss$V1==x),2:ncol(rss)])
  ints = sort.int(y, decreasing = T, index.return = T)
  colnames(rss[,2:ncol(rss)])[ints$ix[1:40]]
}))
mean.activity = lapply(c("OPC-like","astro-like","intermediate","neuro-like"), function(x) {
  cells = df.ss$cell[which(df.ss$ICA_score_states==x)]
  colMeans(TF.AUCg[which(rownames(TF.AUCg) %in% cells),])
})
top.activity = unlist(lapply(seq(length(mean.activity)), function(x) {
  names(sort(mean.activity[[x]], decreasing = T)[1:40])
}))

for(i in 1:2){
  top.corr = rownames(df.corr[[i]][order(df.corr[[i]]$corr_ps),])[c(1:40,(nrow(df.corr[[i]])-39):nrow(df.corr[[i]]))]
  TFs.plot = unique(c(top.rss,top.activity,top.corr))
  
  known.TFs = list(c("NFIA","NFIB","SOX9","EMX2","RNF20"),
                   c("OLIG1","OLIG2","ASCL1","PRRX1","SOX10","NKX6-2"),
                   c("NEUROD","ZEB1","ZEB2","PAX6"))
  known.regulons = lapply(known.TFs, function(x) {
    unlist(lapply(x, function(y) {colnames(TF.AUCg)[grepl(y,colnames(TF.AUCg))]}))
  })
  names(known.regulons) = c("Astrocyte","Oligodendrocyte","Neural")
  TFs.plot = TFs.plot[!(TFs.plot %in% unlist(known.regulons))]
  df.row = data.frame(Lineage = c(rep("None",length(TFs.plot)),
                                  unlist(lapply(seq(length(known.regulons)), 
                                                function(x) {rep(names(known.regulons)[x],length(known.regulons[[x]]))}))
  ))
  TFs.plot = c(TFs.plot, unlist(known.regulons))
  rownames(df.row) = TFs.plot
  cols.list.row = c("#FFFFFF",gg_color_hue(4)[c(1,3,4)])
  names(cols.list.row) = c("None",names(known.regulons))
  row_anno = HeatmapAnnotation(df = df.row, col = list(Lineage = cols.list.row), which = "row")
  
  #add cell state annotation
  cols.list.col = gg_color_hue(4)
  names(cols.list.col) =  c("astro-like","intermediate","OPC-like","neuro-like")
  df.col = data.frame(ICA_score_states = df.ss[sapply(cells,function(x) {which(rownames(df.ss)==x)}),"ICA_score_states"])
  rownames(df.col) = cells
  bot_anno <- HeatmapAnnotation(df = df.col, col = list(ICA_score_states = cols.list.col))
  
  df.toplot.filt = df.toplot[,sapply(TFs.plot, function(x) {which(colnames(df.toplot)==x)})]
  
  col_fun <- colorRamp2(c(0, 0.1, 0.30), c("lightgrey", "orange","blue"))
  ht <- Heatmap(
    t(df.toplot.filt),
    name = "TFactivity",
    col = col_fun,
    cluster_rows = T,
    cluster_columns = F,
    row_gap = unit(0, "mm"),
    column_gap = unit(0, "mm"),
    show_row_names = T,
    row_names_gp = gpar(fontsize = 6),
    show_column_names = F,
    row_title = NULL,
    column_title = "Pseudotime",
    bottom_annotation  = bot_anno,
    right_annotation = row_anno,
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
  pdf(paste0(fig.dir,"/ps",i,"_heatmap_TFactivity_filt_",cancer,".pdf"), width = 15, height = 12, useDingbats = FALSE)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}



#try other method of filtering
#ADD ANYTHING WHERE 10% IS >0.3
top.rss = unlist(lapply(c("OPC-like","astro-like","intermediate","neuro-like"), function(x) {
  y = as.numeric(rss[which(rss$V1==x),2:ncol(rss)])
  ints = sort.int(y, decreasing = T, index.return = T)
  colnames(rss[,2:ncol(rss)])[ints$ix[1:10]]
}))

for(i in 1:2){
  top.corr = rownames(df.corr[[i]][order(df.corr[[i]]$corr_ps),])[c(1:20,(nrow(df.corr[[i]])-19):nrow(df.corr[[i]]))]
  
  TFs.plot = unique(c(top.rss,top.corr))
  
  cols.list.row = c("#FFFFFF",gg_color_hue(4)[c(1,3,4)])
  names(cols.list.row) = c("None",names(known.regulons))
  row_anno = HeatmapAnnotation(df = df.row, col = list(Lineage = cols.list.row), which = "row")
  
  #add cell state annotation
  names(cols.list.col) =  c("astro-like","intermediate","OPC-like","neuro-like")
  df.col = data.frame(ICA_score_states = df.ss[sapply(cells,function(x) {which(rownames(df.ss)==x)}),"ICA_score_states"])
  rownames(df.col) = cells
  bot_anno <- HeatmapAnnotation(df = df.col, col = list(ICA_score_states = cols.list.col))
  
  df.toplot.filt = df.toplot[,sapply(TFs.plot, function(x) {which(colnames(df.toplot)==x)})]
  
  col_fun <- colorRamp2(c(0, 0.1, 0.30), c("lightgrey", "orange","blue"))
  ht <- Heatmap(
    t(df.toplot.filt),
    name = "TFactivity",
    col = col_fun,
    cluster_rows = T,
    cluster_columns = F,
    row_gap = unit(0, "mm"),
    column_gap = unit(0, "mm"),
    show_row_names = T,
    row_names_gp = gpar(fontsize = 6),
    show_column_names = F,
    row_title = NULL,
    column_title = "Pseudotime",
    bottom_annotation  = bot_anno,
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
  pdf(paste0(fig.dir,"/ps",i,"_heatmap_TFactivity_filt2_",cancer,".pdf"), width = 15, height = 12, useDingbats = FALSE)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}










TF.AUCmean.malig = (TF.AUCg.malig+TF.AUCr.malig)/2
#correlation with pseudotime
for(i in 1){
  #compute corelation
  corr_ps = lapply(seq(ncol(TF.AUCmean.malig)), function(x) {
    res = cor.test(df.ss[,which(colnames(df.ss)==paste0("slingPseudotime_",i))],TF.AUCmean.malig[,x], 
                   use = "complete.obs", method = "pearson")
    return(c(res$estimate, res$p.value))
  })
  
  df.corr = do.call("rbind",corr_ps)
  df.corr = as.data.frame(df.corr)
  rownames(df.corr) = colnames(TF.AUCmean.malig)
  colnames(df.corr) = c("corr_ps","pval_1")
  
  Nplot = 6
  for(direc in c("up","down")){
    plot.list = list()
    plot.list2 = list()
    #plot top 5, both expression on trajectory and dot plot of expr v pseudotime
    for(b in 1:Nplot){
      if(direc=="up"){
        gene = rownames(df.corr[order(df.corr$corr_ps,decreasing=T),])[b]
      }else{
        gene = rownames(df.corr[order(df.corr$corr_ps),])[b]
      }
      df3 = data.frame(pseudotime = df.ss[,which(colnames(df.ss)==paste0("slingPseudotime_",i))], 
                       activity = TF.AUCmean.malig[,which(colnames(TF.AUCmean.malig)==gene)],
                       cluster = df.ss$kmeans_DR)
      plot.list2[[b]] = ggplot(df3, aes(x = pseudotime, y = activity, color = cluster)) +
        geom_point(size=pt.sz) + labs(title=gene)
      
      df3 = cbind(df3, DR1=df.ss$DR1, DR2=df.ss$DR2)
      plot.list[[b]] = ggplot(df3, aes(x = DR1, y = DR2, color = activity)) + geom_point(size=pt.sz) #+
      # geom_path(data = df.crv, aes(x=crv1_1,y=crv1_2),col="black") +
      # geom_path(data = df.crv, aes(x=crv2_1,y=crv2_2),col="black")
    }
    p1 = wrap_plots(c(plot.list,plot.list2), ncol = 6, byrow = T)
    print(p1)
    ggsave(paste0(fig.dir,"/ps",i,"_TF_mean_",direc,"_",cancer,"_",j,"_",dr,".pdf"),
           width = 20, height = 5)
  }
}




#positively and negatively correlated TFs separately
TF.AUCdiff.malig = TF.AUCr.malig - TF.AUCg.malig

#start with one regulon
ssrand = sample(seq(nrow(TF.AUCg.malig)),nrow(TF.AUCg.malig))
regs = rownames(df.corr[order(df.corr$corr_ps),])[(nrow(df.corr)-29):nrow(df.corr)]
plot.list = list()
count = 1
for(reg in regs){
  reg = unlist(strsplit(reg,"(", fixed=T))[1]
  ind = which(grepl(reg,colnames(TF.AUCdiff),fixed=T))
  df.diff = data.frame(AUCgene = TF.AUCg.malig[ssrand,ind], AUCregion = TF.AUCr.malig[ssrand,ind], 
                       AUCdiff = TF.AUCdiff.malig[ssrand,ind], cell_state = df.ss$ICA_score_states[ssrand])
  plot.list[[count]] = ggplot(df.diff, aes(x = AUCgene, y = AUCregion, color = cell_state)) + geom_point() +
    ggtitle(reg)
  count = count + 1
}
wrap_plots(plot.list, ncol=5)
ggsave(paste0(fig.dir,"/scatter_AUC_gene_v_region_pos_corr_",cancer,".pdf"), width = 24, height = 24)


regs = rownames(df.corr[order(df.corr$corr_ps),])[1:30]
plot.list = list()
count = 1
for(reg in regs){
  reg = unlist(strsplit(reg,"(", fixed=T))[1]
  ind = which(grepl(reg,colnames(TF.AUCdiff),fixed=T))
  df.diff = data.frame(AUCgene = TF.AUCg.malig[ssrand,ind], AUCregion = TF.AUCr.malig[ssrand,ind], 
                       AUCdiff = TF.AUCdiff.malig[ssrand,ind], cell_state = df.ss$ICA_score_states[ssrand])
  plot.list[[count]] = ggplot(df.diff, aes(x = AUCgene, y = AUCregion, color = cell_state)) + geom_point() +
    ggtitle(reg)
  count = count + 1
}
wrap_plots(plot.list, ncol=5)
ggsave(paste0(fig.dir,"/scatter_AUC_gene_v_region_neg_corr_",cancer,".pdf"), width = 24, height = 24)




###REPEAT PLOTS with mean of gene and reg and diff gene and reg
####Pseudo time plots of well correlated and large difference side by side
#### any difference by grade


##### OLD WAY - TO REMOVE?



# #change column names so they are less long
# colnames(TF.scores) = sapply(colnames(TF.scores), function(x) {
#   paste(unlist(strsplit(x,"_"))[c(1,2,4)],collapse="_")
# })

# reg.spec = read.csv("../../SCENICplus_local/figures/rss_astro_test.csv", sep=";",
#                      header=T, row.names = 1)
# 
# #for each malignant cell type, find top 5 TFs
# malig.cell.types = c("oligo_traj","astro_traj","undifferentiated","neuro_traj")
# inds.to.plot = unlist(lapply(malig.cell.types, function(x) {
#   sorted = sort.int(as.numeric(reg.spec[which(rownames(reg.spec)==x),]),decreasing=T, index.return = T)
#   sorted$ix[1:5]
# }))
# 
# 
# df.to.plot = cbind(df[which(rownames(df) %in% rownames(TF.scores)),], TF.scores[unlist(sapply(rownames(df), function(x) {which(rownames(TF.scores)==x)})),
#                                  inds.to.plot])
# TFs.to.plot = colnames(TF.scores)[inds.to.plot]
# 
# 
# plot.list = list()
# for(i in 1:length(TFs.to.plot)){
#   plot.list[[i]] = ggplot(df.to.plot, aes(x = DR1, y = DR2, color = .data[[TFs.to.plot[i]]])) + geom_point(size=pt.sz) +
#     geom_path(data = df.crv, aes(x=crv1_1,y=crv1_2),col="black") +
#     geom_path(data = df.crv, aes(x=crv2_1,y=crv2_2),col="black")
# }
# p1 = wrap_plots(plot.list, ncol = 5)
# ggsave(paste0(fig.dir,"/top_TFs_",cancer,"_",j,"_",dr,".pdf"), width = 25, height = 15)


#find cells in both dfs
joint.cells = intersect(rownames(df),rownames(TF.scores))
df.ss = df[sapply(joint.cells, function(x) {which(rownames(df)==x)}),]
TF.scores.malig = TF.scores[sapply(joint.cells, function(x) {which(rownames(TF.scores)==x)}),]

#for both pseudotimes
for(i in 1:2){
  #compute corelation
  corr_ps = lapply(seq(ncol(TF.scores.malig)), function(x) {
    res = cor.test(df.ss[,which(colnames(df.ss)==paste0("slingPseudotime_",i))],TF.scores.malig[,x], 
                   use = "complete.obs", method = "pearson")
    return(c(res$estimate, res$p.value))
  })
  
  df.corr = do.call("rbind",corr_ps)
  df.corr = as.data.frame(df.corr)
  rownames(df.corr) = colnames(TF.scores.malig)
  colnames(df.corr) = c("corr_ps","pval_1")
  
  Nplot = 6
  for(direc in c("up","down")){
    plot.list = list()
    plot.list2 = list()
    #plot top 5, both expression on trajectory and dot plot of expr v pseudotime
    for(b in 1:Nplot){
      if(direc=="up"){
        gene = rownames(df.corr[order(df.corr$corr_ps,decreasing=T),])[b]
      }else{
        gene = rownames(df.corr[order(df.corr$corr_ps),])[b]
      }
      df3 = data.frame(pseudotime = df.ss[,which(colnames(df.ss)==paste0("slingPseudotime_",i))], 
                       activity = TF.scores.malig[,which(colnames(TF.scores.malig)==gene)],
                       cluster = df.ss$kmeans)
      plot.list2[[b]] = ggplot(df3, aes(x = pseudotime, y = activity, color = cluster)) +
        geom_point(size=pt.sz) + labs(title=gene)
      
      df3 = cbind(df3, DR1=df.ss$DR1, DR2=df.ss$DR2)
      plot.list[[b]] = ggplot(df3, aes(x = DR1, y = DR2, color = activity)) + geom_point(size=pt.sz) #+
        # geom_path(data = df.crv, aes(x=crv1_1,y=crv1_2),col="black") +
        # geom_path(data = df.crv, aes(x=crv2_1,y=crv2_2),col="black")
    }
    p1 = wrap_plots(c(plot.list,plot.list2), ncol = 6, byrow = T)
    print(p1)
    ggsave(paste0(fig.dir,"/ps",i,"_TF_act_",direc,"_",cancer,"_",j,"_",dr,".pdf"),
           width = 20, height = 5)
  }
}


#find TFs most DE active in each trajectory-based cell state
#plot average activity in pseudotime bins along trajectory
TF.scores.malig$state = sapply(rownames(TF.scores.malig), function(x) {df.ss$cell.state.traj[which(rownames(df.ss)==x)]})

state_means <- df %>%
  group_by(state) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
top_tfs <- state_means %>%
  tidyr::pivot_longer(-state,
                      names_to = "TF",
                      values_to = "activity") %>%
  group_by(state) %>%
  slice_max(activity, n = 10)

#or rescale
state_means_z <- state_means

state_means_z[,-1] <- scale(state_means[,-1])
top_tfs_z <- state_means_z %>%
  pivot_longer(-state) %>%
  group_by(state) %>%
  slice_max(value, n = 10)


#or wilcoxon
library(presto)  # fast Wilcoxon

markers <- wilcoxauc(
  X = as.matrix(df[, -which(colnames(df)=="state")]),
  y = df$state
)
head(markers)


#visualise in heatmap
library(pheatmap)

top <- unique(top_tfs$TF)

mat <- state_means_z[, top]
rownames(mat) <- state_means_z$state

pheatmap(mat, scale="none")


#use tradeseq
library(tradeSeq)

fit <- fitGAM(
  counts = t(df[,-state_col]),
  pseudotime = df$pseudotime,
  cellWeights = df$lineage
)





#plot a few known genes from astro/oligo signature
genes.to.plot = c("DLGAP1","OPCML","FGF12","SOX6","DSCAM","TNR","GRIA4",
                  "DPP10","ADGRV1","LINC01088","AC012405.1","CADPS","SPARCL1")
genes.to.plot = genes.to.plot[genes.to.plot %in% rownames(sce)]

df = df[,1:15]
df = cbind(df, as.data.frame(as.matrix(t(counts(sce)[sapply(genes.to.plot, 
                                                  function(x) {which(rownames(sce)==x)}),]))))

plot.list = list()
for(i in 1:length(genes.to.plot)){
  plot.list[[i]] = ggplot(df, aes(x = DR1, y = DR2, color = .data[[genes.to.plot[i]]])) + geom_point(size=pt.sz) +
    geom_path(data = df.crv, aes(x=crv1_1,y=crv1_2),col="black") +
    geom_path(data = df.crv, aes(x=crv2_1,y=crv2_2),col="black")
}
p1 = wrap_plots(plot.list, ncol = ceiling(length(genes.to.plot)/2))
ggsave(paste0(fig.dir,"/ICA_lin_markers_",cancer,"_",j,"_",dr,".pdf"), width = 18, height = 5)


#markers from Hamed et al
#plot a few known genes
genes.to.plot2 = list(
  NC1 = c("FOXD3","ERBB3","PLP1","SOX10","NGFR","ID3"), #neural crest first half
  NC2 = c("ETS1","CD44","SPARC","HES1","ANXA2","VIM"), #neural crest second half
  MES_SC = c("TNC","MET","ALDH1A3","NES","VGF"), #mes stem cell
  NSC = c("EGFR","SOX9","SLC1A3","ALDOC","NTSR2"), # NSC
  OPC = c("PDGFRA", "OLIG1", "OLIG2", "C1QL1","CSPG4","CNP"), # OPC
  NPC = c("CD24", "SOX11", "SOX4", "DLL3" ,"STMN4"), # NPC
  IMM_OL = c("FYN", "ENPP6", "TNR", "MINK1", "MBP","PLP","APC","MOG","MAG"),  #immature OL
  AC = c("AGT", "GJA1", "AQP4","CLU","APOE") # AC
)
genes.to.plot2 = lapply(genes.to.plot2, function(x) {x[x %in% rownames(sce)]})
df = df[,1:15]
mx_pl = max(sapply(genes.to.plot2, length))

df = cbind(df, as.data.frame(as.matrix(t(counts(sce)[sapply(unlist(genes.to.plot2), 
                                                            function(x) {which(rownames(sce)==x)}),]))))
df2 = df[sample(seq(nrow(df)),nrow(df)),]
plot.list = list()
for(i in 1:length(genes.to.plot2)){
  for(m in 1:mx_pl){
    if(m<=length(genes.to.plot2[[i]])){
      plot.list[[(i-1)*mx_pl + m]] = ggplot(df2, aes(x = DR1, y = DR2, color = .data[[genes.to.plot2[[i]][m]]])) + geom_point(size=pt.sz) +
        geom_path(data = df.crv, aes(x=crv1_1,y=crv1_2),col="black") +
        geom_path(data = df.crv, aes(x=crv2_1,y=crv2_2),col="black")
    }else{
      plot.list[[(i-1)*mx_pl + m]] = plot_spacer()
    }
  }
}
p1 = wrap_plots(plot.list, ncol = mx_pl)
ggsave(paste0(fig.dir,"/cell_type_markers_",cancer,"_",j,"_",dr,".pdf"), width = 18, height = 20)


if(cancer=="astro"){
  #find out what separate cell group is next to clusters 2 and 6 are
  #roughly set a line to separate these cells and define as new cluster
  #use only cluster 6 as this is where bulk of cells are
  inds.to.keep = which(sce$kmeans==6)
  sce2 = sce[,inds.to.keep]
  df2 = df[inds.to.keep,]
  sce2$kmeans[which((0.005 - df2$DR1) > df2$DR2)] = 5
  #find markers
  markers_default2 <- findMarkers(sce2, 
                                  groups = factor(sce2$kmeans), # clusters to compare
                                  test.type = "wilcox",   # t-test (default)
                                  direction = "up", # test for either higher or lower expression (default)
                                  lfc = 0, # null hypothesis log-fold-change = 0 (default)
                                  pval.type = "any", # ranking of p-values based on any comparison (default)
                                  assay.type = "counts"
  )
}




#filter counts to speed up fitGAM
counts_filt = as.matrix(data.sct)[rowSums(as.matrix(data.sct))>1000,]

#fit model
sceGAM <- fitGAM(counts = counts_filt, sds = SlingshotDataSet(sce))

plotGeneCount(sce, clusters = cl2, models = sceGAM)

#test for association with pseudotime for each branch separately
pseudotime_association <- associationTest(sceGAM, global=F, lineages=T)
pseudotime_association <- pseudotime_association[order(pseudotime_association[,"waldStat_1"], decreasing = T), ]
pseudotime_association <- cbind(pseudotime_association, p.adjust(pseudotime_association[,"pvalue_1"], method = "fdr"))
pseudotime_association <- cbind(pseudotime_association, p.adjust(pseudotime_association[,"pvalue_2"], method = "fdr"))
colnames(pseudotime_association) = c(colnames(pseudotime_association)[1:(length(colnames(pseudotime_association))-2)], "fdr_1","fdr_2")
pseudotime_association <- pseudotime_association[order(pseudotime_association[,"pvalue_1"]), ]
pseudotime_association = cbind(pseudotime_association,rownames(pseudotime_association))
colnames(pseudotime_association)[length(colnames(pseudotime_association))] = "feature_id"
feature_id = pseudotime_association[,"feature_id"]

plot_differential_expression <- function(feature_id, n) {
  # feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  plot_row = cowplot::plot_grid(plotGeneCount(sce, gene = feature_id[n], clusters = cl2, models = sceGAM) + 
                       ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sceGAM, counts_filt, gene = feature_id[n]))
  title <- cowplot::ggdraw() +  cowplot::draw_label(feature_id[n],fontface = 'bold',x = 0,hjust = 0) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  cowplot::plot_grid(title, plot_row, ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )
}

plot_gene_expression <- function(geneID) {
  # feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  plot_row = cowplot::plot_grid(plotGeneCount(sce, gene = geneID, clusters = cl2, models = sceGAM) + 
                                  ggplot2::theme(legend.position = "none"), 
                                plotSmoothers(sceGAM, counts_filt, gene = geneID))
  title <- cowplot::ggdraw() +  cowplot::draw_label(geneID, fontface = 'bold', x = 0, hjust = 0) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  cowplot::plot_grid(title, plot_row, ncol = 1,
                     # rel_heights values control vertical title margins
                     rel_heights = c(0.1, 1)
  )
}

plot_differential_expression(feature_id,4)

#association test
pseudotime_association <- associationTest(sceGAM, global=F, lineages=T)
pseudotime_association = as.data.frame(pseudotime_association)
pseudotime_association <- cbind(pseudotime_association, p.adjust(pseudotime_association[,"pvalue_1"], method = "fdr"))
pseudotime_association <- cbind(pseudotime_association, p.adjust(pseudotime_association[,"pvalue_2"], method = "fdr"))
colnames(pseudotime_association) = c(colnames(pseudotime_association)[1:(length(colnames(pseudotime_association))-2)], 
                                     "fdr_1","fdr_2")
#compute correlation between expression and pseudotime
corr_1 = sapply(rownames(pseudotime_association), function(x) {
  cor(slingPseudotime(sce)[,1],counts(sce)[which(rownames(sce)==x),], use = "complete.obs")
})
corr_2 = sapply(rownames(pseudotime_association), function(x) {
  cor(slingPseudotime(sce)[,2],counts(sce)[which(rownames(sce)==x),], use = "complete.obs")
})
pseudotime_association = cbind(pseudotime_association,corr_1,corr_2)

#order by stat
pseudotime_association <- pseudotime_association[order(pseudotime_association$waldStat_1, decreasing = T), ]
feature_id = rownames(pseudotime_association)

featN = 1
if(pseudotime_association$corr_1[featN]>0){
  plot_differential_expression(feature_id,featN)
}


pseudotime_association <- pseudotime_association[order(pseudotime_association$waldStat_2, decreasing = T), ]
feature_id = rownames(pseudotime_association)

featN = 1
if(pseudotime_association$corr_2[featN]>0){
  plot_differential_expression(feature_id,featN)
}


plot_gene_expression("DSCAM")


pseudotime_association <- pseudotime_association[order(pseudotime_association$corr_2, decreasing = T), ]


plotGeneCount(curve = sce, clusters = cl2,title = "Colored by cluster")

pseudotime <- slingPseudotime(crv1, na=FALSE)
cellWeights <- slingCurveWeights(crv1)

icMat2 <- evaluateK(counts = as.matrix(data.sct), pseudotime = pseudotime, cellWeights = cellWeights,
                    k=3:7, nGenes = 100, verbose = FALSE, plot = TRUE)



############
####### old code



# library(Seurat)
# #extract cells which are defined as neuro-like and see where they fall on integrated umap of all samples
# oligo.ss = readRDS(paste0("/Users/alyne/oligo.ss_integrated_RPCA_all_",cancer,".rds"))
# neuro.cells = rownames(df)[which(df$cell.state.traj=="neuro-like")]
# neuro_logi = sapply(names(oligo.ss$orig.ident), function(x) {
#   x %in% neuro.cells
# })
# oligo.ss$neuro_logi = neuro_logi
# 
# dim.red = "umap.fmnn"
# if(dim.red == "umap.fmnn"){
#   dim.str = "fmnn"
#   reduc = "integrated.mnn"
# }else if(dim.red == "umap.har"){
#   dim.str = "rpca"
#   reduc = "integrated.dr"
# }
# #make data frame so can plot SNVs last
# df.SNV = data.frame(umaphar_1 = Embeddings(oligo.ss,dim.red)[,1],
#                 umaphar_2 = Embeddings(oligo.ss,dim.red)[,2], SNV = oligo.ss$IDH.mut, 
#                 cluster = oligo.ss$fmnn_clusters, sample = oligo.ss$orig.ident)
# levs = unique(df.SNV$SNV)[!is.na(unique(df.SNV$SNV))]
# cols = gg_color_hue(length(levs))
# df.SNV$color = sapply(df.SNV$SNV, function(x) {
#   if(is.na(x)){return("gray70")}else{cols[which(levs==x)]}
# })
# #sort to put all NAs first, i.e., to plot all the other things last
# df.SNV = df.SNV[order(df.SNV$SNV, na.last = F),]
# 
# #also add trajectory cell type
# df.SNV$cell.state.traj = sapply(rownames(df.SNV), function(x) {
#   if(x %in% rownames(df)){
#     df$cell.state.traj[which(rownames(df)==x)]
#   }else{
#     NA
#   }
# })
# 
# #also add trajectory cell type
# df.SNV$icnv.score = sapply(rownames(df.SNV), function(x) {
#   if(x %in% rownames(df)){
#     df$icnv.score[which(rownames(df)==x)]
#   }else{
#     NA
#   }
# })
# 
# #add ICA UMAP
# df.umap = do.call("rbind",lapply(rownames(df.SNV), function(x) {
#   if(x %in% rownames(u)){
#     u[which(rownames(u)==x),]
#   }else{
#     c(NA,NA)
#   }
# }))
# colnames(df.umap) = c("UMAP1","UMAP2")
# df.SNV = cbind(df.SNV, df.umap)
# 
# oligo.ss <- FindNeighbors(oligo.ss, reduction = reduc, dims = 1:30, graph.name = paste0(dim.str,"_nn"))
# oligo.ss <- FindClusters(oligo.ss, resolution = 0.3, cluster.name = paste0(dim.str,"_clusters"), graph.name = paste0(dim.str,"_nn"))
# 
# pt.sz = 0.001
# lg.pt.sz = 3
# p1 = DimPlot(oligo.ss, group.by = "cell.type.sr.imm", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
# p2 = DimPlot(oligo.ss, group.by = "iCNV.malig.ident", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
# df.SNV2 = rbind(df.SNV[which(is.na(df.SNV$cell.state.traj)),],df.SNV[which(df.SNV$cell.state.traj!="neuro-like"),],
#                 df.SNV[which(df.SNV$cell.state.traj=="neuro-like"),])
# p3 = ggplot(df.SNV2, aes(x = umaphar_1, y = umaphar_2, color = cell.state.traj)) + geom_point(size=0.1)+
#   guides(colour = guide_legend(override.aes = list(size=lg.pt.sz)))
# p4 = ggplot(df.SNV, aes(x = umaphar_1, y = umaphar_2, color = SNV)) + geom_point(size=0.1)+
#   guides(colour = guide_legend(override.aes = list(size=lg.pt.sz)))
# p5 = DimPlot(oligo.ss, group.by = paste0(dim.str,"_clusters"), reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
# p6 = ggplot(df.SNV, aes(x=cluster, y=icnv.score, fill = cluster)) +
#   geom_boxplot(outlier.shape = NA) + geom_jitter(shape=20, position=position_jitter(0.2),size=0.1)
# p7 = ggplot(df.SNV, aes(x = UMAP1, y = UMAP2, color = cell.state.traj)) + geom_point(size=0.1)+
#   guides(colour = guide_legend(override.aes = list(size=lg.pt.sz)))
# p8 = ggplot(df.SNV, aes(x = UMAP1, y = UMAP2, color = sample)) + geom_point(size=0.1)+
#   guides(colour = guide_legend(override.aes = list(size=lg.pt.sz)))
# p9 = DimPlot(oligo.ss, group.by = "orig.ident", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
# # p10 = FeaturePlot(oligo.ss, features = c("RBFOX3","SLC17A7","GAD2"), reduction = dim.red, pt.size = pt.sz, raster=F, ncol=3)
# p10 = FeaturePlot(oligo.ss, features = c("CLDN5","ADAMTS13","CD93"), reduction = dim.red, pt.size = pt.sz, raster=F, ncol=3)
# print((p1 + p2 + p3 ) / (p4 + p5 + p6) / (p7 + p8 + p9) / p10)
# ggsave(paste0(fig.dir,"/umaps_trajectory_defined_cell_states_",dim.str,".pdf"), height = 14, width = 18)
# 
# cells.to.check = rownames(df.SNV[which(df.SNV$umaphar_1>2 & df.SNV$umaphar_1<4.5 & df.SNV$umaphar_2<2.5 & df.SNV$umaphar_2>-2),])
# oligo.ss <- SetIdent(oligo.ss, cells = cells.to.check, value = 'unknown')
# oligo.ss <- SetIdent(oligo.ss, cells = colnames(oligo.ss)[!(colnames(oligo.ss) %in% cells.to.check)], value = 'everything_else')
# results = FindMarkers(oligo.ss, ident.1 = 'unknown', ident.2 = 'everything_else')
# 
# Idents(oligo.ss) <- paste0(dim.str,"_clusters")
# results = FindMarkers(oligo.ss, ident.1 = 13, ident.2 = 0)
# write.table(rownames(results[results$avg_log2FC>0,]), 
#             file = paste0("unknown_markers_",dim.str,".tsv"), sep="\t", row.names=F, col.names = F)
# 
# #NEFL - neuron marker
# 
# save(df, df.ICA, df.SNV, sce, file = paste0("post_traj_output_",cancer,"_",j,rem,".RData"))
# 
# df.filt = df
# df.ICA.filt = df.ICA
# df.SNV.filt = df.SNV
# 
# load(paste0("pre_traj_output_",cancer,"_",j,".RData"))
# 
# #compute nearest neighbour grph on all samples and return nearest neighbours
# #use rpca integration since fmnn doesn't push all the samples together
# dim.red = "umap.har"
# if(dim.red == "umap.fmnn"){
#   dim.str = "fmnn"
#   reduc = "integrated.mnn"
# }else if(dim.red == "umap.har"){
#   dim.str = "rpca"
#   reduc = "integrated.dr"
# }
# city_neighborhood <- FindNeighbors(oligo.ss, reduction = reduc, dims = 1:30, 
#                           graph.name = paste0(dim.str,"_nn"), return.neighbor=TRUE)
# 
# #add trajectory cells type for all the cells except the removed sample
# #hmm prob need to reduce number of neighbours
# cell.state.traj = character(nrow(df))
# for(i in 1:nrow(df)){
#   x = rownames(df)[i]
#   if(x %in% rownames(df.filt)){
#     cell.state.traj[i] = df.filt$cell.state.traj[which(rownames(df.filt) == x)]
#   }else{
#     ind = which(colnames(oligo.ss)==x)
#     cul_de_sac<-TopNeighbors(city_neighborhood[[paste0(dim.str,"_nn")]],cell=colnames(oligo.ss)[ind],n=20)
#     cs = names(sort(table(df.filt$cell.state.traj[which(rownames(df.filt) %in% cul_de_sac)]), 
#                     decreasing = T, na.last = T))[1]
#     if(is.null(cs)){
#       #find label of nearest point and transfer
#       print(i)
#       print("nearest neighbours all in unlabelled sample")
#       dists = apply(df[which(df$sample!=samps.to.rem),c(1,2)], 1, function(y) {
#         sqrt(sum((y - df[i,1:2])^2))
#       })
#       cell.state.traj[i] = unlist(df$cell.state.traj[which(df$sample!=samps.to.rem)[which.min(dists)]])
#       print(cell.state.traj[i])
#     }else{
#       cell.state.traj[i] = unlist(cs)
#     }
#   }
# }
# 
# df$cell.state.traj = cell.state.traj



## other dimensionality reductions
# #compute PCA on ICA?
# pca_ica = calculatePCA(data.ica,ncomponents = 2)
# reducedDim(sce, "PCA_ICA") <- pca_ica
# 
# #read in LLE dim reduction and add to sce
# i=5
# lle = fread(paste0(cancer,"_samps_ss_rem_cc_",j,"_lle_",i,".tsv.gz"), data.table=F, header=T)
# lle = lle[,2:3]
# colnames(lle) = c("LLE1","LLE2")
# rownames(lle) = colnames(data.sct)
# reducedDim(sce, "LLE") <- lle
#
# #compute phate dim reduction and add to sce object, differernt nn values
# plot.list = list()
# for(i in 2:10){
#   ica.phate = phate(t(data.ica), knn = i)
#   df = data.frame(DR1 = ica.phate$embedding[,1], DR2 = ica.phate$embedding[,2])
#   df = cbind(df, cell.md)
#   plot.list[[i-1]] = ggplot(df, aes(x = DR1, y = DR2, colour = cell.type.sr.imm)) + geom_point()
# }
# wrap_plots(plot.list, ncol=3)
# #compute phate dim reduction on ICA projection and add to sce object
# ica.phate = phate(t(data.ica))
# reducedDim(sce, "PHATE_ICA") <- ica.phate$embedding
# 
# #compute umap
# u = uwot::umap(t(data.ica))
# reducedDim(sce, "UMAP") <- u




## other options of dimension reductions
# #compute PCA on SCT data
# pca = calculatePCA(sce, exprs_values = "counts")
# reducedDim(sce, "PCA") <- pca
# 
# reducedDim(sce, "ICA") <- t(data.ica)
# 
# # other dimensionality reductions
# #read in LLE dim reduction and add to sce
# i=5
# lle = fread(paste0("../astro_ss_rem_cc_1/",
#                    cancer,"_samps_ss_rem_cc_1_lle_",i,".tsv.gz"), data.table=F, header=T)
# lle = lle[,2:3]
# colnames(lle) = c("LLE1","LLE2")
# cell.md.ss = fread(paste0("../astro_ss_rem_cc_1/metadata_",cancer,
#                           "_samps_ss_rem_cc_1.tsv.gz"), data.table=F)
# rownames(cell.md.ss) = cell.md.ss[,1]
# cell.md.ss = cell.md.ss[,-1]
# # rownames(lle) = colnames(data.sct)
# # reducedDim(sce, "LLE") <- lle
# #compute phate dim reduction on ICA projection and add to sce object
# ica.phate = phate(t(data.ica))
# reducedDim(sce, "PHATE_ICA") <- ica.phate$embedding
# #compute umap
# u = uwot::umap(t(data.ica))
# reducedDim(sce, "UMAP") <- u
# df.dr = data.frame(DM1 = dm$DC1, DM2 = dm$DC2,
#                    UMAP1 = u[,1], UMAP2 = u[,2],
#                    PHATE1 = ica.phate$embedding[,1], PHATE2 = ica.phate$embedding[,2])
# df.dr = cbind(df.dr,cell.md)
# df.dr.ss = data.frame(LLE1 = lle[,1], LLE2 = lle[,2])
# df.dr.ss = cbind(df.dr.ss,cell.md.ss)
# #plot figure of different dim reductions
# plot.list = list()
# plot.list[[1]] = ggplot(df.dr, aes(x = DM1, y = DM2, colour = tirosh)) + geom_point(size=pt.sz) + 
#   guides(color = guide_legend(override.aes = list(size=2)))
# plot.list[[2]] = ggplot(df.dr, aes(x = UMAP1, y = UMAP2, colour = tirosh)) + geom_point(size=pt.sz) + 
#   guides(color = guide_legend(override.aes = list(size=2)))
# plot.list[[3]] = ggplot(df.dr, aes(x = PHATE1, y = PHATE2, colour = tirosh)) + geom_point(size=pt.sz) + 
#   guides(color = guide_legend(override.aes = list(size=2)))
# plot.list[[4]] = ggplot(df.dr.ss, aes(x = LLE1, y = LLE2, colour = tirosh)) + geom_point(size=pt.sz) + 
#   guides(color = guide_legend(override.aes = list(size=2)))
# p1 = wrap_plots(plot.list)
# print(p1)
# ggsave(paste0(fig.dir,"/dr_",cancer,"_",j,".pdf"), width = 10, height = 8)










plot.list = list()
#plot LLE with clusters
plot.list[[1]] = ggplot(df, aes(x = DR1, y = DR2, colour = kmeans)) + geom_point(size=pt.sz) + 
  guides(color = guide_legend(override.aes = list(size=2)))

#tirosh cell types
df2 = df[sample(seq(nrow(df)),nrow(df)),]
plot.list[[2]] = ggplot(df2, aes(x = DR1, y = DR2, colour = tirosh.cell.type.new)) + geom_point(size=pt.sz) + 
  guides(color = guide_legend(override.aes = list(size=2)))

#SingleR cell type
#sort df so neurons are at bottom
df2 = rbind(df[!(df$cell.type.sr.imm %in% c("GABAergic","Glutamatergic")),],
            df[df$cell.type.sr.imm %in% c("GABAergic","Glutamatergic"),])
plot.list[[3]] = ggplot(df2, aes(x = DR1, y = DR2, colour = cell.type.sr.imm)) + geom_point(size=pt.sz) + 
  theme(legend.key.size = unit(0.2, 'cm')) + guides(colour = guide_legend(override.aes = list(size=3)))

#proportion cycling for each cluster and cycling cells
prop.cyc.cl = sapply(unique(df$kmeans), function(x) {
  sum(df$cyc.status[which(df$kmeans==x)]=="cycling")/length(df$cyc.status[which(df$kmeans==x)])
})
df$prop.cyc = sapply(df$kmeans, function(x) {
  prop.cyc.cl[which(unique(df$kmeans)==x)]
})
plot.list[[4]] = ggplot(df, aes(x = DR1, y = DR2, colour = prop.cyc)) + geom_point(size=pt.sz)
df2 = rbind(df[!(df$cyc.status %in% c("cycling")),], df[df$cyc.status %in% c("cycling"),])
plot.list[[5]] = ggplot(df2, aes(x = DR1, y = DR2, colour = cyc.status)) + geom_point(size=pt.sz) + 
  guides(color = guide_legend(override.aes = list(size=2)))

#plot cytotrace value
df2 = df[sample(seq(nrow(df)),nrow(df)),]
plot.list[[6]] = ggplot(df2, aes(x = DR1, y = DR2, colour = CytoTRACE)) + geom_point(size=pt.sz)
plot.list[[7]] = ggplot(df2, aes(x = DR1, y = DR2, colour = traj.cell.type.new)) + geom_point(size=pt.sz) + 
  theme(legend.key.size = unit(0.2, 'cm')) + guides(colour = guide_legend(override.aes = list(size=3)))
plot.list[[8]] = ggplot(df2, aes(x = DR1, y = DR2, colour = CytoTRACE2)) + geom_point(size=pt.sz)

#sample and ini/rec and icnv score
plot.list[[9]] = ggplot(df2, aes(x = DR1, y = DR2, colour = CytoTRACE2.pot)) + geom_point(size=pt.sz)
plot.list[[10]] = ggplot(df2, aes(x = DR1, y = DR2, colour = sample)) + geom_point(size=pt.sz) + 
  theme(legend.key.size = unit(0.2, 'cm')) + guides(colour = guide_legend(override.aes = list(size=3)))
plot.list[[11]] = ggplot(df2, aes(x = DR1, y = DR2, colour = ini_rec)) + geom_point(size=pt.sz) + 
  guides(color = guide_legend(override.aes = list(size=2)))
plot.list[[12]] = ggplot(df2, aes(x = DR1, y = DR2, colour = icnv.score)) + geom_point(size=pt.sz)
plot.list[[13]] = ggplot(df2, aes(x = DR1, y = DR2, colour = as.factor(grade))) + geom_point(size=pt.sz) + 
  guides(color = guide_legend(override.aes = list(size=2)))
plot.list[[14]] = ggplot(df2, aes(x = DR1, y = DR2, colour = treatment)) + geom_point(size=pt.sz) + 
  guides(color = guide_legend(override.aes = list(size=2)))
plot.list[[15]] = ggplot(df2, aes(x = DR1, y = DR2, colour = atlas.label)) + geom_point(size=pt.sz) + 
  theme(legend.key.size = unit(0.2, 'cm')) + guides(colour = guide_legend(override.aes = list(size=3)))
plot.list[[16]] = ggplot(df2, aes(x = DR1, y = DR2, colour = gest_week)) + geom_point(size=pt.sz) + 
  guides(color = guide_legend(override.aes = list(size=2)))

p1 = wrap_plots(plot.list, ncol=4)
print(p1)
ggsave(paste0(fig.dir,"/pre_traj_",cancer,"_",j,"_",dr,rem,".pdf"), width = 18, height = 12)







# colour.by = df$cell.type.sr.imm
# cols.uni = gg_color_hue(length(unique(colour.by)))
# cols.tir = cols.uni[factor(colour.by)]
# plot3d(x=dm$DC1, y=dm$DC2, z=dm$DC3, col = cols.tir)
# legend3d("topright", legend = levels(factor(colour.by)), pch = 16, col = cols.uni, cex=1.5)
# 
# open3d()
# colour.by = df$grade
# cols.uni = gg_color_hue(length(unique(colour.by)))
# cols.tir = cols.uni[factor(colour.by)]
# plot3d(x=dm$DC1, y=dm$DC2, z=dm$DC3, col = cols.tir)
# legend3d("topright", legend = levels(factor(colour.by)), pch = 16, col = cols.uni, cex=1.5)
# 
# 
# cl3 <- kmeans(cbind(dm$DC1,dm$DC2,dm$DC3), centers = 7)$cluster
# df$kmeans3 = cl3
# 
# open3d()
# colour.by = df$kmeans3
# cols.uni = gg_color_hue(length(unique(colour.by)))
# cols.tir = cols.uni[factor(colour.by)]
# plot3d(x=dm$DC1, y=dm$DC2, z=dm$DC3, col = cols.tir)
# legend3d("topright", legend = levels(factor(colour.by)), pch = 16, col = cols.uni, cex=1.5)