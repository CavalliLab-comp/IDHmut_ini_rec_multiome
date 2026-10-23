library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(readxl)
library(tidyr)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

setwd("/mnt/beegfs/home/alyne/stagein")
fig.dir = "/mnt/beegfs/home/alyne/stageout/figures_cohort_summary"

samp.md <- as.data.frame(read_excel("../persistent/sample_summary_astro_oligo_anon.xlsx", col_types = "text"))


#read in seurat objects with all samples
astro = readRDS("scenicplus_astro/oligo.ss_integrated_RPCA_all_astro_doub_rem.rds")
oligo = readRDS("scenicplus_oligo/oligo.ss_integrated_RPCA_all_oligo_doub_rem.rds")


#set cell type of malignant cells
astro$cell.type.sr.cl.imm.postmerge[which(astro$iCNV.malig.ident=="malignant")] = "Malignant"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$iCNV.malig.ident=="malignant")] = "Malignant"

#rename endo-related cell types to Endothelial
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="VLMC")] = "Endo"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="VLMC")] = "Endo"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Endo")] = "Endothelial"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Endo")] = "Endothelial"

#rename
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Astro")] = "Astrocyte"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Astro")] = "Astrocyte"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Oligo")] = "Oligodendrocyte"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Oligo")] = "Oligodendrocyte"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Macrophage")] = "Myeloid"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Macrophage")] = "Myeloid"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="GABAergic")] = "GABAergic neuron"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="GABAergic")] = "GABAergic neuron"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Glutamatergic")] = "Glutamatergic neuron"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Glutamatergic")] = "Glutamatergic neuron"

#rename lymphoid
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge %in% 
                                            c("T-cells","B-cells","NK cells"))] = "Lymphoid"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge %in% 
                                            c("T-cells","B-cells","NK cells"))] = "Lymphoid"


#plot full UMAPs of each sample integrated via three different approaches, version in paper is fmnn
sz=0.001
p1 = DimPlot(astro, group.by = "cell.type.sr.cl.imm.postmerge", reduction = "umap.unintegrated",
                 pt.size = sz, raster=F) + theme(text = element_text(size = 20)) +
  labs(x="UMAP1",y="UMAP2",title="Astrocytoma") + NoLegend()
p2 = DimPlot(oligo, group.by = "cell.type.sr.cl.imm.postmerge", reduction = "umap.unintegrated",
                 pt.size = sz, raster=F) + theme(text = element_text(size = 20)) +
  labs(x="UMAP1",y="UMAP2",title="Oligodendroglioma")
print(p1 + p2)
ggsave(paste0(fig.dir,"/umaps_full_both_cancers_unintegrated.pdf"), width = 12, height = 5)

p1 = DimPlot(astro, group.by = "cell.type.sr.cl.imm.postmerge", reduction = "umap.rpca",
             pt.size = sz, raster=F) + theme(text = element_text(size = 20)) +
  labs(x="UMAP1",y="UMAP2",title="Astrocytoma") + NoLegend()
p2 = DimPlot(oligo, group.by = "cell.type.sr.cl.imm.postmerge", reduction = "umap.rpca",
             pt.size = sz, raster=F) + theme(text = element_text(size = 20)) +
  labs(x="UMAP1",y="UMAP2",title="Oligodendroglioma")
print(p1 + p2)
ggsave(paste0(fig.dir,"/umaps_full_both_cancers_rpca.pdf"), width = 12, height = 5)

p1 = DimPlot(astro, group.by = "cell.type.sr.cl.imm.postmerge", reduction = "umap.har",
             pt.size = sz, raster=F) + theme(text = element_text(size = 20)) +
  labs(x="UMAP1",y="UMAP2",title="Astrocytoma") + NoLegend()
p2 = DimPlot(oligo, group.by = "cell.type.sr.cl.imm.postmerge", reduction = "umap.har",
             pt.size = sz, raster=F) + theme(text = element_text(size = 20)) +
  labs(x="UMAP1",y="UMAP2",title="Oligodendroglioma")
print(p1 + p2)
ggsave(paste0(fig.dir,"/umaps_full_both_cancers_har.pdf"), width = 12, height = 5)

p1 = DimPlot(astro, group.by = "cell.type.sr.cl.imm.postmerge", reduction = "umap.fmnn",
             pt.size = sz, raster=F) + theme(text = element_text(size = 20)) +
  labs(x="UMAP1",y="UMAP2",title="Astrocytoma") + NoLegend()
p2 = DimPlot(oligo, group.by = "cell.type.sr.cl.imm.postmerge", reduction = "umap.fmnn",
             pt.size = sz, raster=F) + theme(text = element_text(size = 20)) +
  labs(x="UMAP1",y="UMAP2",title="Oligodendroglioma")
print(p1 + p2)
ggsave(paste0(fig.dir,"/umaps_full_both_cancers_fmnn.pdf"), width = 12, height = 5)

DimPlot(astro, group.by = "cell.type.sr.cl.imm.postmerge", reduction = "umap.fmnn",
             pt.size = sz, raster=F) + theme(text = element_text(size = 20),axis.text.x=element_blank(),axis.text.y=element_blank()) +
  labs(x="UMAP1",y="UMAP2",title="Astrocytoma") + NoLegend()
ggsave(paste0(fig.dir,"/umaps_full_astro_fmnn.pdf"), width = 5, height = 5)

DimPlot(oligo, group.by = "cell.type.sr.cl.imm.postmerge", reduction = "umap.fmnn",
             pt.size = sz, raster=F) + theme(text = element_text(size = 20),axis.text.x=element_blank(),axis.text.y=element_blank()) +
  labs(x="UMAP1",y="UMAP2",title="Oligodendroglioma")
ggsave(paste0(fig.dir,"/umaps_full_oligo_fmnn.pdf"), width = 8, height = 5)


#############
## make barplot summarising cell type proportions
#############

df.astro = data.frame(table(astro$cell.type.sr.cl.imm.postmerge), cancer = "Astrocytoma")
df.oligo = data.frame(table(oligo$cell.type.sr.cl.imm.postmerge), cancer = "Oligodendroglioma")
df.cell.types = rbind(df.astro,df.oligo)

cols = gg_color_hue(length(unique(df.cell.types$Var1)))
names(cols) = sort(unique(df.cell.types$Var1))

df.cell.types$Var1 = factor(df.cell.types$Var1, levels = df.astro$Var1[order(df.astro$Freq)])

ggplot(df.cell.types, aes(fill=Var1, y=Freq, x=cancer)) + 
  geom_bar(position="fill", stat="identity")+ labs(fill='Cell type',x = "", y = "Proportion") +
  theme(text=element_text(size=20)) + scale_fill_manual(values = cols)
ggsave(paste0(fig.dir,"/barplot_celltypes_subtype.pdf"))


