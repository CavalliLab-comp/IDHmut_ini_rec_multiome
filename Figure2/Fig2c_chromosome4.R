library(Seurat)
library(Signac)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


cancer = "oligo"
data.dir = paste0("/mnt/beegfs/home/alyne/stagein/scenicplus_",cancer)


#read in seurat object with malignant cells from all smples within subtype
fig.str = "malig"
if(cancer=="astro"){
  oligo.ss = readRDS(paste0(data.dir,"/oligo.ss_integrated_RPCA_",fig.str,"_",cancer,"_doub_rem.rds"))
}else{
  oligo.ss = readRDS(paste0(data.dir,"/oligo.ss_integrated_RPCA_",fig.str,"_",cancer,".rds"))
}


###oligo sample which clearly have chr4 deletion in all cells
#oligo_P0003_rec1_S0006
#oligo_P0022_ini_S0042

###clonal
#oligo_P0007_ini_S0015
#oligo_P0024_rec1_S0047


#get sample id
oligo.ss$sample_id = sapply(oligo.ss$orig.ident, function(x) {
  unlist(strsplit(x,"_"))[2]
})


chr4_del = c("oligo_P0003_rec1_S0006","oligo_P0022_ini_S0042")
chr4_no_del = c("oligo_P0001_ini_S0001","oligo_P0001_rec1_S0002","oligo_P0001_rec2_S0003","oligo_P0002_ini_S0004",
               "oligo_P0002_rec1_S0005","oligo_P0004_ini_S0009","oligo_P0004_rec1_S0010","oligo_P0005_ini_S0011",
               "oligo_P0005_rec1_S0012","oligo_P0008_ini_S0017","oligo_P0010_ini_S0021","oligo_P0012_ini_S0023",
               "oligo_P0012_rec1_S0024","oligo_P0024_ini_S0046","oligo_P0030_rec1_S0059","oligo_P0031_ini_S0060",
               "oligo_P0031_rec1_S0061")

#keep sample with either clera do or do not have chr 4 deletion
samps.to.keep = sapply(c(chr4_del,chr4_no_del), function(x) {
  unlist(strsplit(x,"_"))[4]
})
chr4_del_id = sapply(chr4_del, function(x) {
  unlist(strsplit(x,"_"))[4]
})
chr4_no_del_id = sapply(chr4_no_del, function(x) {
  unlist(strsplit(x,"_"))[4]
})

oligo.ss.filt = subset(oligo.ss, subset = sample_id %in% samps.to.keep)

#make dataframe of sample id, ICA17 gene set module score, and chr4 deletion status
df = oligo.ss.filt@meta.data[,c("sample_id","ICA_17_1")]
df$chr4 = sapply(df$sample_id, function(x) {
  if(x %in% chr4_del_id){
    "del"
  }else{
    "no_del"
  }
})
df$chr4 = factor(df$chr4)

#plot histogram
ggplot(df, aes(x=ICA_17_1, fill=chr4, color=ICA_17_1)) + 
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins = 100) +
  labs(x = "chr4 metacomponent score") + theme(text = element_text(size = 16)) + lims(x = c(-0.2,0.4))
ggsave("histogram_chr4_score_by_group.pdf", height=6, width = 5)


