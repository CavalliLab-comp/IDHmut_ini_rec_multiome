#script to make figures summarising cell types and changes ini/rec
#cell type proportions are transformed using chiPower to allow use of standard
#statistical tests on compositional data (Greenacre 2024)

library(dplyr)
library(patchwork)
library(ggplot2)
library(Seurat)
library(Signac)
library(readxl)
library(tidyr)
library(tibble)
library(ggpubr)
library(rstatix)
library(textshape)
library(compositions)
library(easyCODA)
library(vegan)

#function to get ggplot colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

setwd("/mnt/beegfs/home/alyne/stagein")
fig.dir = '/mnt/beegfs/home/alyne/stageout/figures_cohort_summary'
dir.create(fig.dir)

# #calcsub version
# setwd("/data/users/alyne/gliomult/svn/analyse/results")
# fig.dir = '/data/users/alyne/gliomult/svn/analyse/results/figures_cohort_summary'
# dir.create(fig.dir)

#read in full seurat objects with all samples for both subtypes
astro = readRDS("scenicplus_astro/oligo.ss_integrated_RPCA_all_astro_doub_rem.rds")
oligo = readRDS("scenicplus_oligo/oligo.ss_integrated_RPCA_all_oligo_doub_rem.rds")

#add patient as meta data
astro$patient = sapply(colnames(astro), function(x) {unlist(strsplit(x,"_"))[2]})
oligo$patient = sapply(colnames(oligo), function(x) {unlist(strsplit(x,"_"))[2]})
astro$sample = sapply(colnames(astro), function(x) {unlist(strsplit(x,"_"))[4]})
oligo$sample = sapply(colnames(oligo), function(x) {unlist(strsplit(x,"_"))[4]})



#rename cell types to broad categories/cleaner names
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Endo")] = "Endothelial"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="VLMC")] = "Endothelial"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Astro")] = "Astrocyte"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Oligo")] = "Oligodendrocyte"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="T-cells")] = "Lymphoid"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Macrophage")] = "Myeloid"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="GABAergic")] = "GABAergic neuron"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Glutamatergic")] = "Glutamatergic neuron"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Endo")] = "Endothelial"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="VLMC")] = "Endothelial"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Astro")] = "Astrocyte"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Oligo")] = "Oligodendrocyte"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="T-cells")] = "Lymphoid"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="B-cells")] = "Lymphoid"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="NK cells")] = "Lymphoid"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Macrophage")] = "Myeloid"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="GABAergic")] = "GABAergic neuron"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Glutamatergic")] = "Glutamatergic neuron"

#rename malignant
astro$cell.type.sr.cl.imm.postmerge[which(astro$iCNV.malig.ident=="malignant")] = "Malignant"
df.astro = data.frame(table(astro$cell.type.sr.cl.imm.postmerge), cancer = "Astrocytoma")


#rename malignant
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$iCNV.malig.ident=="malignant")] = "Malignant"
df.oligo = data.frame(table(oligo$cell.type.sr.cl.imm.postmerge), cancer = "Oligodendroglioma")



#make combined data frame
df.cell.types = rbind(df.astro,df.oligo)
df.cell.types$cancer = factor(df.cell.types$cancer)

#sort by order of abundance in astro
df.cell.types$Var1 = factor(df.cell.types$Var1, levels = df.astro$Var1[order(df.astro$Freq)])


#make colour palette
cellt.pal = gg_color_hue(length(unique(df.cell.types$Var1)))
names(cellt.pal) = sort(as.character(unique(df.cell.types$Var1)))


#make paired figures showing changes ini/rec
#keep non-abundant cell types but we won't plot them

df.astro.pat = data.frame(cell.type = astro$cell.type.sr.cl.imm.postmerge, ini_rec = astro$ini_rec,
                          patient = astro$patient, cancer = "Astrocytoma", sample = astro$sample)
#compute proportions in each sample
df_counts.astro <- df.astro.pat %>%
  count(patient, ini_rec, cancer, cell.type, sample, name = "count") %>%
  group_by(sample) %>%
  complete(cell.type = unique(df.astro.pat$cell.type),
    fill = list(count = 0)) %>%
  fill(patient, ini_rec, cancer, .direction = "downup") %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()


#same for oligo
df.oligo.pat = data.frame(cell.type = oligo$cell.type.sr.cl.imm.postmerge, ini_rec = oligo$ini_rec,
                          patient = oligo$patient, cancer = "Oligodendroglioma", sample = oligo$sample)
df_counts.oligo <- df.oligo.pat %>%
  count(patient, ini_rec, cancer, cell.type, sample, name = "count") %>%
  group_by(sample) %>%
  complete(cell.type = unique(df.oligo.pat$cell.type),
    fill = list(count = 0)) %>%
  fill(patient, ini_rec, cancer, .direction = "downup") %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

df_counts = rbind(df_counts.astro,df_counts.oligo)
df_counts$cancer = factor(df_counts$cancer)


#########
## chiPower transformation
#########

#make count matrix of samples x cell type counts
gli.counts0 <- df_counts %>%
  select(sample, cell.type, count) %>%
  pivot_wider(names_from = cell.type, values_from = count, values_fill = 0) %>%
  column_to_rownames("sample") %>%
  as.matrix()
gli.counts = gli.counts0 +1

### --------------------------------------------------------------------------------
### chiPower function
chiPower<- function(data, close=TRUE, power=1, chi=TRUE, BoxCox=TRUE, center=TRUE) {
  # function to compute chiPower transformation with various options
  # close: close data
  # power: tranformation (i.e., value of lambda)
  # chi: chi-square standardization
  # BoxCox: apply 1/power rescaling (without the subtraction of 1)
  # center: center the final result by column means
  foo <- data^power
  if(close)   foo <- foo / rowSums(foo)
  if(chi)     foo <- sweep(foo, 2, sqrt(colMeans(foo)), FUN="/")
  if(BoxCox)  foo <- (1/power) * foo
  if(center)  foo <- sweep(foo, 2, colMeans(foo))
  if(!center) foo <- foo - 1/sqrt(ncol(data))
  data.chiPower <- foo
  data.chiPower
}
### ----------------------------------------------------------------------------------

### logratio distances (with column weight constant in this case, i.e. weight=FALSE)
### Use CLRs to get logratio distances
gli.clr <- CLR(CLOSE(gli.counts), weight=FALSE)$LR
### LRA is PCA of CLRs (or you can get it directly using function LRA)
gli.lra <- PCA(gli.clr, weight=FALSE)
### Principal row coordinates define logratio geometry (last dimension has no variance)
gli.lra.rpc <- gli.lra$rowpcoord[,-(ncol(gli.counts)-1)]
### Logratio distances from CLRs
gli.clr.dist <- dist(gli.clr)/sqrt(ncol(gli.counts))


### loop on powers from 0.001 to 1 and find optimal using Procrustes
### comparison with chiPower transforms in two different ways
### 1. Spearman rank correlation between distances
### 2. Procrustes correlation between geometries
npower <- 100
gli.cor <- rep(0,npower)
gli.pro <- rep(0,npower)  
### Original composition on data set with 0s
gli0.comp <- CLOSE(gli.counts0)    
### start of loop (takes a while) -----------
for(power in 1:npower) {
  foo.chiPower <- chiPower(gli0.comp, power=power/npower)
  foo.dist     <- dist(foo.chiPower)
  foo.rpc      <- PCA(foo.chiPower, weight=FALSE)$rowpcoord[,1:(ncol(gli0.comp)-1)]
  gli.cor[power] <- cor(gli.clr.dist, foo.dist, method="spearman")
  gli.pro[power] <- protest(gli.lra.rpc, foo.rpc, permutations=0)$t0
}

### plot the Procrustes (black) and Spearman correlations (blue)
pdf(paste0(fig.dir,"/chiPower_correlations_all_malig_combined.pdf"))
plot(gli.pro, type="l", lwd=2, xaxt="n", bty="n", ylim=c(0.3,1),
     xlab="Power of chiPower transformation", ylab="Procrustes correlation / Spearman")
axis(1, at=seq(0,100,10), labels=seq(0, 1, 0.1))
### blue curve for Spearman correlation
lines(gli.cor, col="blue")
dev.off()


#transform counts using Procrustes max
gli.chiPower <- chiPower(gli0.comp, power=max(gli.pro))

#add transformed values to df_counts
df_counts <- df_counts %>% mutate(chiPower = gli.chiPower[cbind(sample, cell.type)])


#relabel rec2 from P0001 to plot as ini/rec
#this makes a second pair from patient P0001, relabelled to P0001b, of rec1 and rec2
#NB this means P0001 rec1 is plotted twice, once as a rec and once as an ini

df_counts.S002 = df_counts[df_counts$sample=="S0002",]
df_counts.S002$ini_rec = "ini"
df_counts.S002$sample = paste0(df_counts.S002$sample,"b")
df_counts.S002$patient = paste0(df_counts.S002$patient,"b")

df_counts$sample[which(df_counts$sample=="S0003")] = paste0(df_counts$sample[which(df_counts$sample=="S0003")],"b")
df_counts$patient[which(df_counts$sample=="S0003b")] = paste0(df_counts$patient[which(df_counts$sample=="S0003b")],"b")
df_counts$ini_rec[which(df_counts$sample=="S0003b")] = "rec1"
df_counts = rbind(df_counts,df_counts.S002)

#add column so can label this pair as different
df_counts$comparison = "ini_rec"
df_counts$comparison[which(grepl("b",df_counts$sample))] = "rec1_rec2"


#keep only paired for ini-rec plots
df_filtered <- df_counts %>%
  group_by(patient) %>%
  filter(all(c("ini", "rec1") %in% ini_rec)) %>%
  ungroup()


#add treatment groups to each patient
samp.md = as.data.frame(read_excel("~/persistent/sample_summary_astro_oligo_anon.xlsx"))
df_filtered$treatment = sapply(df_filtered$patient, function(x) {
  if(x=="P0001b"){
    rel.rows = samp.md[which(samp.md$patient=="P0001"),]
    #get rec1 row and add treatment group
    tr = rel.rows[which(rel.rows$ini_rec=="rec2"),c("TMZ","PCV","Radiotherapy")]
    paste(names(tr[which(tr==1)]), collapse="_")
  }else{
    rel.rows = samp.md[which(samp.md$patient==x),]
    #get rec1 row and add treatment group
    tr = rel.rows[which(rel.rows$ini_rec=="rec1"),c("TMZ","PCV","Radiotherapy")]
    paste(names(tr[which(tr==1)]), collapse="_")
  }
})




#### ini-rec plots for each subtype by treatment group
#### plot proportions, not chiPower
#### not enough patients in each group for statistical significance

ct.to.plot = c("Astrocyte","GABAergic neuron","Glutamatergic neuron","Oligodendrocyte","Myeloid","Malignant")
cellt.pal.filt = cellt.pal[ct.to.plot]
#paired initial/recurrence plots by treatment groups
plot.list2 = list()
for(i in 1:length(ct.to.plot)){
  df.plot <- df_filtered %>%
    filter(cell.type == rev(ct.to.plot)[i]) %>%
    filter(cancer == "Astrocytoma") %>%
    arrange(patient,ini_rec)
  df.plot$ini_rec = factor(df.plot$ini_rec, levels=c("ini","rec1"))
  df.plot$patient = factor(df.plot$patient)
  
  plot.list2[[i]] = ggplot(df.plot, aes(x = ini_rec, y = prop)) + geom_boxplot(outlier.shape=NA,fill = rev(cellt.pal.filt)[i]) +
    geom_line(aes(group = patient)) + facet_wrap(~ treatment) +
    geom_point(size = 2) + labs(y="Proportion",x="") + theme(text = element_text(size = 20)) +
    ggtitle(paste0(rev(ct.to.plot)[i]))
}
p1 = wrap_plots(plot.list2)
print(p1)
ggsave(paste0(fig.dir,"/cell_type_proportions_ini_rec_paired_by_treatment_astro.pdf"), width = 15, height = 10)


plot.list2 = list()
for(i in 1:length(ct.to.plot)){
  df.plot <- df_filtered %>%
    filter(cell.type == rev(ct.to.plot)[i]) %>%
    filter(cancer == "Oligodendroglioma") %>%
    arrange(patient,ini_rec)
  df.plot$ini_rec = factor(df.plot$ini_rec, levels=c("ini","rec1"))
  df.plot$patient = factor(df.plot$patient)
  
  plot.list2[[i]] = ggplot(df.plot, aes(x = ini_rec, y = prop)) + 
    geom_boxplot(outlier.shape=NA,fill = rev(cellt.pal.filt)[i]) +
    geom_line(aes(group = patient, color = comparison)) + facet_wrap(~ treatment) +
    scale_colour_manual(values=c(ini_rec="black",rec1_rec2="chartreuse4")) +
    geom_point(size = 2) + labs(y="Proportion",x="") + theme(text = element_text(size = 20)) +
    ggtitle(paste0(rev(ct.to.plot)[i]))
}
p1 = wrap_plots(plot.list2)
print(p1)
ggsave(paste0(fig.dir,"/cell_type_proportions_ini_rec_paired_by_treatment_oligo.pdf"), width = 20, height = 15)



#refine/combine treatment groups for paper plots


#refine treatment groups
#for astrocytoma, combine all patients with chemo or radio into one group, with second group for untreated
#plot proportions, not enough data in each group for statistical analysis

#paired initial/recurrence plots by treatment groups
plot.list2 = list()
for(i in 1:length(ct.to.plot)){
  df.plot <- df_filtered %>%
    filter(cell.type == rev(ct.to.plot)[i]) %>%
    filter(cancer == "Astrocytoma") %>% filter(treatment != "TMZ") %>%
    arrange(patient,ini_rec)
  df.plot$treatment[df.plot$treatment==""] = "None"
  df.plot$treatment_comb = sapply(df.plot$treatment, function(x) {
    if(grepl("Radiotherapy",x)){
      "Chemo/RT"
    }else{x}
  })
  df.plot$treatment_comb = factor(df.plot$treatment_comb, levels=c("None","Chemo/RT"))
  df.plot$ini_rec = factor(df.plot$ini_rec, levels=c("ini","rec1"))
  df.plot$patient = factor(df.plot$patient)
  
  # paired Wilcoxon per facet
  stat.test <- df.plot %>% group_by(treatment_comb) %>%
    wilcox_test(chiPower ~ ini_rec, paired = TRUE) %>%
    add_significance("p") %>% add_xy_position(x = "ini_rec")
  
  plot.list2[[i]] = ggplot(df.plot, aes(x = ini_rec, y = prop)) + geom_boxplot(outlier.shape=NA, fill = rev(cellt.pal.filt)[i]) +
    geom_line(aes(group = patient, color = treatment),linewidth=1) + facet_wrap(~ treatment_comb, ncol=2) +
    geom_point(size = 1.2, aes(color = treatment)) + labs(y="Proportion",x="") + theme(text = element_text(size = 18)) +
    ggtitle(paste0(rev(ct.to.plot)[i])) + scale_colour_manual(values=c(None="black",PCV_Radiotherapy="darkorange2",TMZ_Radiotherapy="royalblue3"))
  if(i %% 3 != 0){plot.list2[[i]] = plot.list2[[i]] + guides(color="none")}
}
p1 = wrap_plots(plot.list2)
print(p1)
ggsave(paste0(fig.dir,"/cell_type_proportions_ini_rec_paired_by_treatment_comb_astro.pdf"), width = 15, height = 8)



#refine treatment groups for oligo
#put patients into one of three groups: no treatment between ini/rec, PCV or TMZ either with or without radio
#statistical tests on chiPower, but plot proportions

plot.list2 = list()
for(i in 1:length(ct.to.plot)){
  df.plot <- df_filtered %>% filter(cell.type == rev(ct.to.plot)[i]) %>%
    filter(cancer == "Oligodendroglioma") %>% arrange(patient,ini_rec)
  df.plot$treatment[df.plot$treatment==""] = "None"
  df.plot$treatment_comb = sapply(df.plot$treatment, function(x) {
    if(grepl("Radiotherapy",x)){
      "TMZ/RT"
    }else{
      if(grepl("PCV|TMZ",x)){
        "Chemo"
      }else{
        x
      }
    }
  })
  df.plot$treatment_comb = factor(df.plot$treatment_comb, levels = c("None","Chemo","TMZ/RT"))
  df.plot$ini_rec = factor(df.plot$ini_rec, levels=c("ini","rec1"))
  df.plot$patient = factor(df.plot$patient)
  df.plot$comparison[which(df.plot$treatment=="PCV" & df.plot$comparison=="ini_rec")] = "PCV"
  df.plot$comparison[which(df.plot$treatment=="PCV" & df.plot$comparison=="rec1_rec2")] = "PCV rec1_rec2"
  df.plot$comparison[which(df.plot$treatment=="TMZ")] = "TMZ"
  df.plot$chemotherapy = df.plot$comparison
  # paired Wilcoxon per facet
  stat.test <- df.plot %>% group_by(treatment_comb) %>% wilcox_test(chiPower ~ ini_rec, paired = TRUE) %>%
    add_significance("p") %>% add_xy_position(x = "ini_rec")
  
  plot.list2[[i]] = ggplot(df.plot, aes(x = ini_rec, y = prop)) + 
    geom_boxplot(outlier.shape=NA,fill = rev(cellt.pal.filt)[i]) +
    geom_line(aes(group = patient, color = chemotherapy),linewidth=1) + facet_wrap(~ treatment_comb) +
    scale_colour_manual(values=c(ini_rec="black",PCV="chocolate3","PCV rec1_rec2"="orange",TMZ="royalblue3"), breaks = c("PCV","PCV rec1_rec2","TMZ")) +
    geom_point(size = 1.2, aes(color = chemotherapy)) + labs(y="Proportion",x="") + theme(text = element_text(size = 18)) +
    ggtitle(paste0(rev(ct.to.plot)[i]))
  if(i %% 3 != 0){plot.list2[[i]] = plot.list2[[i]] + guides(color="none")}
}
p1 = wrap_plots(plot.list2)
print(p1)
ggsave(paste0(fig.dir,"/cell_type_proportions_ini_rec_paired_by_treatment_comb_oligo.pdf"), width = 15, height = 8)






#supp figures of ini rec proportions for each subtype ignoring treatment groups
#plot proportions, compute statistical significance using wilcoxon on chiPower


#IDH-A
plot.list2 = list()
for(i in 1:length(ct.to.plot)){
  df.plot <- df_filtered %>%
    filter(cell.type == rev(ct.to.plot)[i]) %>%
    filter(cancer == "Astrocytoma") %>%
    arrange(patient,ini_rec)
  df.plot$ini_rec = factor(df.plot$ini_rec, levels=c("ini","rec1"))
  df.plot$patient = factor(df.plot$patient)
  #Wilcoxon test
  pval = wilcox.test(df.plot$chiPower[seq(1,nrow(df.plot),2)], df.plot$chiPower[seq(2,nrow(df.plot),2)], paired = TRUE)$p.value
  
  plot.list2[[i]] = ggplot(df.plot, aes(x = ini_rec, y = prop)) + geom_boxplot(outlier.shape=NA,fill = rev(cellt.pal.filt)[i]) +
    geom_line(aes(group = patient),linewidth=1) +
    geom_point(size = 1.2) + labs(y="Proportion",x="") + theme(text = element_text(size = 14),plot.title = element_text(size = 12)) +
    ggtitle(paste0(rev(ct.to.plot)[i]," (p=", signif(pval, digits=2),")"))
}
p1 = wrap_plots(plot.list2)
print(p1)
ggsave(paste0(fig.dir,"/cell_type_proportions_ini_rec_paired_astro.pdf"), width = 8, height = 6)


#IDH-O
plot.list2 = list()
for(i in 1:length(ct.to.plot)){
  df.plot <- df_filtered %>%
    filter(cell.type == rev(ct.to.plot)[i]) %>%
    filter(cancer == "Oligodendroglioma") %>%
    arrange(patient,ini_rec)
  df.plot$ini_rec = factor(df.plot$ini_rec, levels=c("ini","rec1"))
  df.plot$patient = factor(df.plot$patient)
  #Wilcoxon test
  pval = wilcox.test(df.plot$chiPower[seq(1,nrow(df.plot),2)], df.plot$chiPower[seq(2,nrow(df.plot),2)], paired = TRUE)$p.value
  
  plot.list2[[i]] = ggplot(df.plot, aes(x = ini_rec, y = prop)) + geom_boxplot(outlier.shape=NA,fill = rev(cellt.pal.filt)[i]) +
    geom_line(aes(group = patient, color = comparison),linewidth=1) +
    scale_colour_manual(values=c(ini_rec="black",rec1_rec2="chartreuse4")) +
    geom_point(size = 1.2) + labs(y="Proportion",x="") + theme(text = element_text(size = 14),plot.title = element_text(size = 12)) +
    ggtitle(paste0(rev(ct.to.plot)[i]," (p=", signif(pval, digits=2),")"))
    if(i %% 3 != 0){plot.list2[[i]] = plot.list2[[i]] + guides(color="none")}
}
p1 = wrap_plots(plot.list2)
print(p1)
ggsave(paste0(fig.dir,"/cell_type_proportions_ini_rec_paired_oligo.pdf"), width = 10, height = 6)





#finally, do comparison of proportions in astrocytoma and oligodendroglioma
#remove rec1/rec2 pair from P0001

#remove S0002b and correct S0003b and P0001b
df_counts = df_counts[which(df_counts$sample!="S0002b"),]
df_counts$sample[which(df_counts$sample=="S0003b")] = "S0003"
df_counts$patient[which(df_counts$patient=="P0001b")] = "P0001"
df_counts$ini_rec[which(df_counts$sample=="S0003")] = "rec2"


#figures to compare astro oligo
plot.list2 = list()
for(i in 1:length(ct.to.plot)){
  df.plot <- df_counts %>%
    filter(cell.type == rev(ct.to.plot)[i])
  df.plot$cancer = factor(df.plot$cancer, levels=c("Astrocytoma","Oligodendroglioma"))
  #Wilcoxon test
  pval = wilcox.test(df.plot$chiPower[df.plot$cancer=="Astrocytoma"], df.plot$chiPower[df.plot$cancer=="Oligodendroglioma"])$p.value
  
  plot.list2[[i]] = ggplot(df.plot, aes(x = cancer, y = prop)) + geom_boxplot(outlier.shape=NA,fill = rev(cellt.pal.filt)[i]) +
    labs(y="Proportion",x="") + theme(text = element_text(size = 20)) + geom_jitter(width = 0.15) +
    ggtitle(paste0(rev(ct.to.plot)[i]," (p=", signif(pval, digits=2),")"))
}
p1 = wrap_plots(plot.list2)
print(p1)
ggsave(paste0(fig.dir,"/cell_type_proportions_astro_oligo.pdf"), width = 15, height = 10)



