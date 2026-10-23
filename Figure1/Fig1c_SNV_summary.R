#script to make figures summarising cell types and changes ini/rec
library(dplyr)
library(patchwork)
library(ggplot2)
library(Seurat)
library(Signac)
library(readxl)
library(tidyr)
library(data.table)
library(ComplexHeatmap)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# setwd("/mnt/beegfs/home/alyne/stagein")
# fig.dir = '/mnt/beegfs/home/alyne/stageout/figures_cohort_summary'
# dir.create(fig.dir)
data.dir = "/data/users/alyne/gliomult/svn/analyse/results"
fig.dir = paste0(data.dir,'/figures_cohort_summary')
dir.create(fig.dir)

#read in sample list
samps = sort(read.table("~/scrnaseq/oligo_10x_multiome/oligo_astro_samples_with_malig.txt", header=F)[,2])

#read in meta data
samp.md <- as.data.frame(read_excel(paste0(data.dir,"/sample_summary_astro_oligo_anon.xlsx"), col_types = "text"))


#read in output from cb_sniffer of SNVs detected in each sample
snp.info = list()
for(i in 1:length(samps)){
  snp.info[[i]] = read.table(paste0(data.dir, "/analysis_", samps[i], "/cb_sniffer_", samps[i],
                                    "/", samps[i], "_RNA_counts_CB.tsv"), sep = "\t", header = T,
                             colClasses = c("character",rep("numeric",2),
                                            rep("character",5),rep("numeric",3)))
  if(file.exists(paste0(data.dir, "/analysis_", samps[i], "/cb_sniffer_", samps[i],
                        "/", samps[i], "_ATAC_counts_CB.tsv"))){
    snp.info[[i]] = rbind(snp.info[[i]], read.table(paste0(data.dir, "/analysis_", samps[i], "/cb_sniffer_", samps[i],
                                                           "/", samps[i], "_ATAC_counts_CB.tsv"), sep = "\t", header = T,
                                                    colClasses = c("character",rep("numeric",2),
                                                                   rep("character",5),rep("numeric",3))))
  }
}

#for each sample, extract mutated genes
mut.uni = lapply(snp.info, function(x) {
  unique(x[,c(1,2,7)])
})

#get the indices of mutations in mut.uni
indMut = lapply(seq(length(mut.uni)), function(y) {
  z = mut.uni[[y]]
  cb.sn = snp.info[[y]]
  unlist(sapply(seq(nrow(z)), function(x) {
    rel.rows = cb.sn[cb.sn$chrm==z[x,1] & cb.sn$start==z[x,2],]
    if(sum(rel.rows$alt_count)>0){
      return(x)
    }else{
      return(NULL)
    }
  }))
})

#get mutation list
mut.list = lapply(seq(length(mut.uni)), function(x) {
  mut.uni[[x]][indMut[[x]],3]
})
names(mut.list) = samps
all.muts = unique(unlist(mut.list))
#sort according to abundance
all.muts.sorted = sort(table(unlist(mut.list)),decreasing=T)


mut.mat = do.call("cbind",lapply(mut.list, function(x) {
  case_when(names(all.muts.sorted) %in% x ~ 1,
            !(names(all.muts.sorted) %in% x) ~ 0)
}))
rownames(mut.mat) = names(all.muts.sorted)

#split columns by patient
split_groups = sapply(colnames(mut.mat), function(x) {unlist(strsplit(x,"_"))[2]})
split_groups = factor(split_groups, levels = unique(split_groups))


#annotate columns according to cancer, ini/rec and grade
col_anno_df = samp.md[sapply(colnames(mut.mat), function(x) {which(samp.md$Sample==x)}),c("ini_rec","Grade","cancer","patient")]
rownames(col_anno_df) = colnames(mut.mat)
colnames(col_anno_df)[2:3] = c("grade","subtype")
col_anno_df$subtype[which(col_anno_df$subtype=="astro")] = "Astrocytoma"
col_anno_df$subtype[which(col_anno_df$subtype=="oligo")] = "Oligodendroglioma"

#make colour scales
cols.grade = c("2" = "skyblue", "3" = "royalblue", "4" = "darkblue")
cols.pat = gg_color_hue(length(unique(col_anno_df$patient)))
#order pats by cancer and ini_rec
pat.ord = unique(col_anno_df$patient[order(col_anno_df$subtype, col_anno_df$patient, col_anno_df$ini_rec)])
names(cols.pat) = c(matrix(pat.ord,ncol=2,byrow=T))
cols.ini_rec = c("ini" = "violet", "rec1" = "purple", "rec2" = "purple4")
cols.cancer = c("Astrocytoma" = "darkgreen", "Oligodendroglioma" = "lightgreen")
full.color.scale = list(cols.ini_rec, cols.grade, cols.cancer, cols.pat)
names(full.color.scale) = colnames(col_anno_df)
top_anno <- HeatmapAnnotation(df = col_anno_df, col = full.color.scale, show_legend = FALSE)
col_fun = c("1" = "grey30", "0" = "gray95")


#make heatmap showing whether mutation detected or not in each sample

ht <- Heatmap(
  as.matrix(mut.mat),
  col = col_fun,
  name = "SNV",
  cluster_rows = F,
  cluster_columns = F,
  column_split = col_anno_df[,c("subtype","patient")],
  # column_gap = unit(0.4, "mm"),
  bottom_annotation  = top_anno,
  show_row_names = T,
  show_column_names = T,
  row_title = NULL,
  column_title = NULL
)

lgd_text <- gpar(fontsize = 14)
lgd_title <- gpar(fontsize = 15)

#make all the legends
lgd.ini_rec <- Legend(title = "ini_rec", at = names(cols.ini_rec), legend_gp = gpar(fill = cols.ini_rec), 
                      labels_gp = lgd_text, title_gp = lgd_title)
lgd.grade <- Legend(title = "grade",at = names(cols.grade), legend_gp = gpar(fill = cols.grade),
                    labels_gp = lgd_text, title_gp = lgd_title)
lgd.subtype <- Legend(title = "subtype", at = names(cols.cancer), legend_gp = gpar(fill = cols.cancer),
                      labels_gp = lgd_text, title_gp = lgd_title)
lgd.patient <- Legend(title = "patient",at = names(cols.pat),
                      legend_gp = gpar(fill = cols.pat), labels_gp = lgd_text, title_gp = lgd_title)
lgd.snv <- Legend(title = "SNV", at = c("1", "0"),legend_gp = gpar(fill = c("grey30", "gray95")),
                  labels_gp = lgd_text,title_gp = lgd_title)

# -----------------------------
# Pack annotation legends
# -----------------------------
# max_height forces them into multiple columns

annotation_lgd <- packLegend(lgd.ini_rec, lgd.grade, lgd.subtype,
                             lgd.patient, direction = "vertical", max_height = unit(10, "cm"), column_gap = unit(1, "cm"))
all_lgd <- packLegend(lgd.snv, annotation_lgd,
                      direction = "horizontal", gap = unit(2, "cm"))

pdf(paste0(fig.dir, "/SNV_summary_full_cohort.pdf"), width = 20, height = 7, useDingbats = FALSE )
grid.newpage()
pushViewport(
  viewport(
    layout = grid.layout(
      nrow = 1,
      ncol = 2,
      widths = c(0.75, 0.25)
    )
  )
)
# draw heatmap
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht,show_heatmap_legend = FALSE,
     show_annotation_legend = FALSE, newpage = FALSE)
popViewport()

# add legends
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(all_lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center", "center"))
popViewport()
popViewport()
dev.off()

