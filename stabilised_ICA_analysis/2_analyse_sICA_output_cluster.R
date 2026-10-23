#file to correlate ICA components from each sample within subtype

rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
data.dir = as.character(args[1])
samp.file = as.character(args[2])
script.dir = as.character(args[3])
cancer = as.character(args[4])

#script to analyse sICA output
library(data.table)
library(gplots)
library(ggplot2)
library(cowplot)
library(patchwork)
library(magrittr)
library(dplyr)
library(Seurat)
library(Signac)
library(igraph)
library(fgsea)
library(org.Hs.eg.db)
library(reactome.db)
library(clusterProfiler)
library(msigdbr)
library(ggpubr)
library(mstknnclust)
library(UCell)

setwd(data.dir)

source(paste0(script.dir,"/corrected_GSEA_function.R"))

#read in full list of samples
sample.info = read.table(samp.file, header = F, colClasses = c("character","character","numeric","numeric"), 
                         stringsAsFactors = F, fill = T)
samps = sample.info[,2]
samps = sort(samps)
#remove samples with no malignant cells
samps = samps[!(samps %in% c("oligo_P0003_rec2_S0007","astro_P0014_ini_S0026","oligo_P0022_rec1_S0043","oligo_P0030_ini_S0058"))]


#sample and correlation parameters to use depending on subtype
if(cancer!="all"){
  samps = samps[grepl(cancer,samps)]
  if(cancer=="astro"){
    cor.thresh = 0.45
    minIC = 3
  }else if(cancer=="oligo"){
    cor.thresh = 0.5
    minIC = 5
  }
}

#can filter to take only cells which definitely have CNA
fig.dir = paste0("figures_sICA_",cancer,"_",cor.thresh)
dir.create(fig.dir)

#function to make ggplot colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#Ture Flase palette
cols.TF = gg_color_hue(2)


#for each sample, read in seurta object and sICA output: metasamples
oligo = list()
metasamples = list()
for(i in 1:length(samps)){
  print(i)
  if(file.exists(paste0("analysis_",samps[i],"/oligo_malig_",samps[i],".rds"))){
    oligo[[i]] = readRDS(paste0("analysis_",samps[i],"/oligo_malig_",samps[i],".rds"))
    DefaultAssay(oligo[[i]]) = "SCT"
    if("ATAC" %in% names(oligo[[i]]@assays)){
      oligo[[i]][['ATAC']] = NULL
      oligo[[i]][["lsi"]] <- NULL
      oligo[[i]][["umap.atac"]] <- NULL
      oligo[[i]][["wnn.umap"]] <- NULL
    }
    oligo[[i]][["RNA"]] <- NULL

    #load scores for each metagene for each cell
    metasamples[[i]] = fread(paste0("analysis_",samps[i],"/sICA_",samps[i],"/metasamples_",samps[i],".tsv.gz"),
                             sep = "\t", header = T, data.table = F)
    rownames(metasamples[[i]]) = metasamples[[i]][,1]
    metasamples[[i]] = metasamples[[i]][,2:ncol(metasamples[[i]])]

    cnames = paste0("mg_",seq(0,(ncol(metasamples[[i]])-1)),
                    "_",samps[i])
    colnames(metasamples[[i]]) = sapply(cnames, function(y) {
      paste(unlist(strsplit(y,"_",fixed = T))[c(1,2,6)],collapse="_")
    })
  }else{
    metasamples[[i]] = NA
    oligo[[i]] = NA
  }
}

#remove NAs (shouldn't be any...)
samps = samps[which(!is.na(oligo))]
oligo = oligo[which(!is.na(oligo))]
metasamples = metasamples[which(!is.na(metasamples))]

print(length(samps))


# read in gene scores for each sample
scores = list()
for(i in 1:length(samps)){
  scores[[i]] = fread(paste0("analysis_",samps[i],"/sICA_",samps[i],"/metagenes_",samps[i],".tsv.gz"),
                      data.table = F, header = T, colClasses = "numeric", sep = "\t")
  rownames(scores[[i]]) = scores[[i]][,1]
  scores[[i]] = scores[[i]][,2:ncol(scores[[i]])]
}

#get genes for each sample
genes.list = lapply(scores, function(x) {rownames(x)})

#get intersection of all genes
genes.intersect = Reduce(intersect, genes.list)

#extract only genes in intersect
scores.intersect = lapply(scores, function(x) {
  x = x[sapply(rownames(x), function(y) {y %in% genes.intersect}),]
  #reorder x so that rows are same order as gene list
  x[match(genes.intersect,rownames(x)),]
})

#make one matrix with all gene scores from all components for each sample
scores.intersect = lapply(seq(length(scores.intersect)), function(x) {
  cnames = paste0("mg_",seq(0,(ncol(scores.intersect[[x]])-1)),"_",samps[x])
  colnames(scores.intersect[[x]]) = sapply(cnames, function(y) {
    paste(unlist(strsplit(y,"_",fixed = T))[c(1,2,6)],collapse="_")
  })
  return(scores.intersect[[x]])
})
scores.int.all = do.call(cbind,scores.intersect)


#correlate geps to find which are observed across multiple samples
#compute correlation
corr.meth = "pearson"

#use all genes
corr = cor(scores.int.all, method = corr.meth)

plot_corr = corr
to_keep = apply(plot_corr, 1, function(x) {sum(x>cor.thresh)>1})
plot_corr = plot_corr[to_keep,to_keep]
png(paste0(fig.dir,"/heatmap_",corr.meth,"_correlation_all_components_",cancer,"_filt.png"),
    width = 2000, height = 2000)
heatmap.2(plot_corr, scale = "none", trace = "none",
          cexRow = 1, cexCol = 1, col = heat.colors(110),
          breaks=c(seq(-1,-0.1,length.out=10),seq(0,1,length.out=101)), margins = c(7,7))
dev.off()



#plot linked subgraphs which pass correlation threshold
stat.ao = sapply(colnames(scores.int.all), function(x) {
  samp = unlist(strsplit(x,"_",fixed=T))[3]
  unlist(strsplit(samps[grepl(samp,samps,fixed=T)],"_",fixed=T))[1]
})
stat.pat = sapply(colnames(scores.int.all), function(x) {
  samp = unlist(strsplit(x,"_",fixed=T))[3]
  unlist(strsplit(samps[grepl(samp,samps,fixed=T)],"_",fixed=T))[2]
})
stat.ir = sapply(colnames(scores.int.all), function(x) {
  samp = unlist(strsplit(x,"_",fixed=T))[3]
  unlist(strsplit(samps[grepl(samp,samps,fixed=T)],"_",fixed=T))[3]
})
cor.str = paste(unlist(strsplit(as.character(cor.thresh),".",fixed = T)),collapse="_")
cor.plot = corr
diag(cor.plot) = 0 #as don't want self-loops
to_keep = apply(cor.plot, 1, function(x) {sum(x>cor.thresh)>0})
cor.plot = cor.plot[to_keep,to_keep]
stat.ao = stat.ao[to_keep]
stat.ir = stat.ir[to_keep]
stat.pat = stat.pat[to_keep]
cor.plot = apply(cor.plot, 1, function(x) {
  x[which(x<=cor.thresh)]=0
  return(x) })
ig = graph_from_adjacency_matrix(cor.plot, mode="undirected", weighted=TRUE)


#get all the connected graphs
#use correlation directly as seems to give more connected graph
sub.graphs = table(components(ig)$membership)
sub.graphs.con = names(sub.graphs[sub.graphs>=minIC])
sub.graph.mem = lapply(sub.graphs.con, function(x) {
  names(components(ig)$membership[components(ig)$membership==x])
})

#save a data.frame
mx.ln = max(sapply(sub.graph.mem, length))
sub.graph.mem2 = lapply(sub.graph.mem, function(x) {c(x, rep(NA,mx.ln-length(x)))})
df.mem = as.data.frame(sub.graph.mem2, col.names = paste0("ICA_",1:length(sub.graphs.con)))
write.table(df.mem, file = paste0(fig.dir,"/sub_graph_members_",cancer,".tsv"),
            sep = "\t", quote = F, row.names = F)




#make graph plot using final threshold coloured by metacomponent
cor.str = paste(unlist(strsplit(as.character(cor.thresh),".",fixed = T)),collapse="_")
stat.sg = sapply(colnames(scores.int.all), function(x) {
  colnames(df.mem)[which(apply(df.mem,2,function(y) {x %in% y}))]
})
col.pal = gg_color_hue(ncol(df.mem))
cor.plot = corr
diag(cor.plot) = 0 #as don't want self-loops
to_keep = apply(cor.plot, 1, function(x) {sum(x>cor.thresh)>0})
cor.plot = cor.plot[to_keep,to_keep]
stat.sg = stat.sg[to_keep]
cor.plot = apply(cor.plot, 1, function(x) {
  x[which(x<=cor.thresh)]=0
  return(x) })
ig = graph_from_adjacency_matrix(cor.plot, mode="undirected", weighted=TRUE)
E(ig)$width <- 10*E(ig)$weight
V(ig)$label.cex = 1.2
V(ig)$stat.sg = sapply(stat.sg, function(x) {x[1]})
subgraphs = colnames(df.mem)
V(ig)$color <- sapply(V(ig)$stat.sg, function(x) {
  if(x %in% subgraphs){
    col.pal[which(subgraphs==x)]
  }else{
    NA
  }
})
png(paste0(fig.dir,"/",corr.meth,"_corr_network_",cor.str,"_corr_thr_sg_",cancer,".png"), width = 2000, height = 2000)
plot(ig, vertex.size=5)
legend('topright', legend = subgraphs, col = col.pal, pch=19, cex=2)
dev.off()





#in astro, split ICA_7 (OPC-like) into three groups as this is a massive metacomponent with clear subgroups
#can avoid averaging across different types of OPC-like state
h = heatmap.2(plot_corr, scale = "none", trace = "none", cexRow = 1, cexCol = 1, col = heat.colors(110),
              breaks=c(seq(-1,-0.1,length.out=10),seq(0,1,length.out=101)), margins = c(7,7))
hclust.obj = as.hclust(h$colDendrogram)

if(cancer=="oligo"){
  #no splitting required
}else{
  #find members
  sg = "ICA_7"
  Nsplit = 3
  sg_members = df.mem[[sg]][!is.na(df.mem[[sg]])]
  Ngroups=1
  i=0
  while(Ngroups<Nsplit){
    i=i+1
    groups = cutree(hclust.obj,i)
    sg_groups = sapply(sg_members, function(x) {groups[which(names(groups)==x)]}, USE.NAMES=F)
    Ngroups = length(unique(sg_groups))
  }
  for(j in 1:Nsplit){
    to.add = names(sg_groups)[which(sg_groups==unique(sg_groups)[j])]
    df.mem[[paste0(sg,"_",j)]] = c(to.add, rep(NA,nrow(df.mem)-length(to.add)))
  }
  df.mem[[sg]] = NULL
}
#rewrite subgraph members
write.table(df.mem, file = paste0(fig.dir,"/sub_graph_members_split_",cancer,".tsv"),
            sep = "\t", quote = F, row.names = F)




#make final plot with split subgroups
cor.str = paste(unlist(strsplit(as.character(cor.thresh),".",fixed = T)),collapse="_")
stat.sg = sapply(colnames(scores.int.all), function(x) {
  colnames(df.mem)[which(apply(df.mem,2,function(y) {x %in% y}))]
})
col.pal = gg_color_hue(ncol(df.mem))
cor.plot = corr
diag(cor.plot) = 0 #as don't want self-loops
to_keep = apply(cor.plot, 1, function(x) {sum(x>cor.thresh)>0})
cor.plot = cor.plot[to_keep,to_keep]
stat.sg = stat.sg[to_keep]
cor.plot = apply(cor.plot, 1, function(x) {
  x[which(x<=cor.thresh)]=0
  return(x) })
ig = graph_from_adjacency_matrix(cor.plot, mode="undirected", weighted=TRUE)
E(ig)$width <- 10*E(ig)$weight
V(ig)$label.cex = 1.2
V(ig)$stat.sg = sapply(stat.sg, function(x) {x[1]})
subgraphs = colnames(df.mem)
V(ig)$color <- sapply(V(ig)$stat.sg, function(x) {
  if(x %in% subgraphs){
    col.pal[which(subgraphs==x)]
  }else{
    NA
  }
})
png(paste0(fig.dir,"/",corr.meth,"_corr_network_",cor.str,
           "_corr_thr_sg_",cancer,"_split.png"), width = 2000, height = 2000)
plot(ig, vertex.size=5)
legend('topright', legend = subgraphs, col = col.pal, pch=19, cex=2)
dev.off()






#for each connected component, find average score for each gene,
#then do gsea and enrichment analysis

m_t2g = msigdbr(species = "Homo sapiens") %>% dplyr::select(gs_name, gene_symbol)

#which MSigDB databases to use
cats = c("H","C8","C2")


#for each subgraph define an 'average' metacomponent
#keep distribution the same
av = "med"
mg.av = list()
for(i in 1:ncol(df.mem)){
  #select all members of this metacomponent
  members = df.mem[,i][which(!is.na(df.mem[,i]))]
  #get relevant column from score matrix
  cols = sapply(members, function(x) {which(colnames(scores.int.all)==x)})
  #define 'average'
  #find median value for each rank
  scores.int.all.st = apply(scores.int.all[,cols], 2, function(y) {
    sort(y, decreasing = T)
  })
  av.distr = rowMedians(scores.int.all.st)
  #find median rank for each gene
  ranks = apply(scores.int.all[,cols], 2, function(x) {
    length(x) - rank(x) + 1
  })
  ranks.av = rank(rowMedians(ranks), ties.method = "random")

  #match gene names based on rank to each value in 'average' metcomponent
  mg.av[[i]] = av.distr[ranks.av]
  names(mg.av[[i]]) = rownames(scores.int.all)
}
df.mg.av = as.data.frame(mg.av)
colnames(df.mg.av) = colnames(df.mem)
fwrite(df.mg.av, file = paste0(fig.dir,"/all_ICA_components_",av,"_",cancer,".tsv.gz"),
       sep = "\t", row.names = T, col.names = T)

#write top and bottom 100 genes from each component
top.genes.med = apply(df.mg.av, 2, function(x) {names(sort(x, decreasing = T))[1:100]})
bott.genes.med = apply(df.mg.av, 2, function(x) {names(sort(x, decreasing = F))[1:100]})
write.table(top.genes.med, paste0(fig.dir,"/top_gene_",av,"_",cancer,".tsv"),
            row.names = F, col.names = T, quote = F, sep="\t")
write.table(bott.genes.med, paste0(fig.dir,"/bott_gene_",av,"_",cancer,".tsv"),
            row.names = F, col.names = T, quote = F, sep="\t")



#get matrix for projecting into ICA space for trajectory analysis
#each column is 'average' component from above (needs to compute inverse for projection)
#NB all components kept including cell cycle
mixing.mat = do.call(cbind, mg.av)
fwrite(as.data.frame(mixing.mat),
       file = paste0(fig.dir,"/mixing_mat_rem_cc_",cancer,".tsv.gz"),
       sep = "\t", append = F, col.names = T, quote = F)


# run GSEA and enrichment analysis using GO and MSigDB
# good to use multiple databases as some components are lineage-based and some are function-based

eGOBP.high.list = list()
eMSIGDB.high.list = list()
#iterate through each metacomponent
for(i in 1:ncol(df.mem)){
  #extract individual members
  sub.graph = colnames(df.mem)[i]
  members = df.mem[,i][which(!is.na(df.mem[,i]))]
  #get relevant column from score matrix
  cols = sapply(members, function(x) {which(colnames(scores.int.all)==x)})
  
  #compute median rank for each gene
  ranks = apply(scores.int.all[,cols], 2, function(x) {length(x) - rank(x) + 1})
  rank.meds = apply(ranks, 1, median)
  
  #compute enrichment in top 100 genes, keep thresholds low to get broad picture
  nUse = 100
  pval.thr = 0.1
  qval.thr = 0.1
  eGOBP.high.list[[i]] = enrichGO(names(sort(rank.meds, decreasing = F)[1:nUse]), OrgDb ="org.Hs.eg.db",
                                  keyType ="SYMBOL", ont ="BP", pvalueCutoff = pval.thr, minGSSize = 10, 
                                  maxGSSize = 500, qvalueCutoff = qval.thr, universe = genes.intersect)
  if(sum(eGOBP.high.list[[i]]@result$pvalue<pval.thr &
         eGOBP.high.list[[i]]@result$p.adjust<pval.thr &
         eGOBP.high.list[[i]]@result$qvalue<qval.thr)>0){
    png(paste0(fig.dir,"/sICA_enrichment_GOBP_top_",nUse,
               "_",sub.graph,"_",cancer,".png"), height = 750, width = 650)
    print(barplot(eGOBP.high.list[[i]], showCategory = 8, label_format = 50, font.size = 20,title = sub.graph))
    dev.off()
  }

  eMSIGDB.high.list[[i]] = enricher(names(sort(rank.meds, decreasing = F)[1:nUse]), TERM2GENE = m_t2g, universe = genes.intersect,
                                    pvalueCutoff = pval.thr, minGSSize = 10, maxGSSize = 500, qvalueCutoff = qval.thr)
  if(sum(eMSIGDB.high.list[[i]]@result$pvalue<pval.thr &
         eMSIGDB.high.list[[i]]@result$p.adjust<pval.thr &
         eMSIGDB.high.list[[i]]@result$qvalue<qval.thr)>0){
    png(paste0(fig.dir,"/sICA_enrichment_MSIGDB_top_",nUse,
               "_",sub.graph,"_",cancer,".png"), height = 750, width = 650)
    print(barplot(eMSIGDB.high.list[[i]], showCategory = 20, label_format = 25, font.size = 16, title = sub.graph))
    dev.off()
  }

  #plot distributions of scores for cells in each sample
  metasamples.rel = lapply(metasamples, function(x) {
    x[,colnames(x) %in% names(cols)]
  })
  metasamples.rel = lapply(names(cols), function(x) {
    ind = which(sapply(metasamples, function(y) {
      x %in% colnames(y)
    }))
    ind2 = which(colnames(metasamples[[ind]])==x)
    return(metasamples[[ind]][,ind2])
  })
  names(metasamples.rel) = names(cols)
  plot.list = lapply(seq(length(metasamples.rel)), function(x) {
    qplot(metasamples.rel[[x]], geom="histogram") + ggtitle(names(metasamples.rel)[x]) +
      theme(axis.title.x=element_blank())
  })
  p1 = wrap_plots(plot.list)
  print(p1)
  ggsave(paste0(fig.dir,"/hists_cell_scores_",sub.graph,"_",cancer,".pdf"), width = 15, height = 15)
  
  #GSEA
  #get median rank for each gene 
  ranks = apply(scores.int.all[,cols], 2, function(x) {
    length(x) - rank(x) + 1
  })
  rank.meds = apply(ranks, 1, median)
  rank.ints = rank(rank.meds, ties.method = "random")
  #use one sample's scores as input to gsea
  #choose the one with largest mean correlation to all other components
  comp.to.use = which.max(colSums(corr[cols,cols]))
  rank.comp = nrow(scores.int.all) - rank(scores.int.all[,comp.to.use],ties.method = "random") + 1
  gsea.input = sapply(seq(nrow(scores.int.all)), function(x) {
    scores.int.all[which(rank.comp==rank.ints[x]),comp.to.use]
  })
  names(gsea.input) = rownames(scores.int.all)

  barplot.list = list()
  count = 1
  for(cat in cats){
    if(cat=="C2"){
      msigdbr_df <- msigdbr(species = "human", category = cat, subcategory = "CP:WIKIPATHWAYS")
    }else{
      msigdbr_df <- msigdbr(species = "human", category = cat)
    }
    pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
    #remove repeated genes within pathways as this was causing a problem in GSEA
    pathwaysH = lapply(pathwaysH, unique)

    fgseaRes = fgseaCellIDcorr(pathways=pathwaysH, gsea.input, minSize = 15, maxSize = 500)

    #tidy results
    fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))

    #Plot the normalized enrichment scores.
    #Color the bar indicating whether or not the pathway was significant
    #show top and bottom n
    fgseaResPos = fgseaResTidy %>% filter(NES > 0)
    fgseaResNeg = fgseaResTidy %>% filter(NES < 0)
    nPrint = 10
    if(nrow(fgseaResPos)>=(nPrint)){
      fgseaResPos = fgseaResPos[1:nPrint,]
    }
    if(nrow(fgseaResNeg)>=(nPrint)){
      fgseaResNeg = fgseaResNeg[(nrow(fgseaResNeg)-nPrint+1):nrow(fgseaResNeg),]
    }
    fgseaResTidy = bind_rows(fgseaResPos, fgseaResNeg)
    fgseaResTidy$psign = sign(fgseaResTidy$NES)*(-log10(fgseaResTidy$padj))
    fgseaResTidy$significance = fgseaResTidy$padj<0.05
    fgseaResTidy$TFcolor = sapply(fgseaResTidy$significance, function(x) {if(x){cols.TF[2]}else{cols.TF[1]}})
    
    barplot.list[[count]] = ggplot(fgseaResTidy, aes(reorder(pathway, psign), NES, fill = I(TFcolor))) +
      geom_col() + coord_flip() + theme_minimal() + labs(x="Pathway", y="Normalized Enrichment Score", title=cat)
    count = count + 1
  }
  p1 = wrap_plots(barplot.list, ncol = 2, nrow = 2)
  print(p1)
  ggsave(paste0(fig.dir,"/sICA_GSEA_MSIGDB_",sub.graph,"_",cancer,".pdf"),
         width = 17, height = 10)
}






###################
## now score all cells to define score-based malignant states
###################






