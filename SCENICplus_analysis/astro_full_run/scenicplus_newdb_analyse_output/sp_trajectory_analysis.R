library(data.table)
library(ggplot2)
library(patchwork)

fig.dir = "/Users/alyne/Library/CloudStorage/Dropbox/Curie/Curie_StCloud/anne-marie/10x_RNA_ATAC_oligo_analysis/SCENICplus_soroban/astro_full_run/figures_scenicplus_astro"

cancer = "astro"
setwd(paste0("/Users/alyne/Library/CloudStorage/Dropbox/Curie/Curie_StCloud/anne-marie/",
             "10x_RNA_ATAC_oligo_analysis/trajectory_analysis_local/",cancer,"_ss_rem_cc_all"))

sp.dir = paste0("/Users/alyne/Library/CloudStorage/Dropbox/Curie/Curie_StCloud/anne-marie",
                "/10x_RNA_ATAC_oligo_analysis/SCENICplus_soroban/",cancer,"_full_run/outs")

#read in TF activities
TF_act.dir = fread(paste0(sp.dir,"/direct_gene_based_AUC_",cancer,".tsv.gz"), data.table = F)
TF_act.ex = fread(paste0(sp.dir,"/extended_gene_based_AUC_",cancer,".tsv.gz"), data.table = F)
TF_act = cbind(TF_act.dir, TF_act.ex[,2:ncol(TF_act.ex)])
#change cell names to match trajectory
TF_act$Cell = sapply(TF_act$Cell, function(x) {unlist(strsplit(x,"___"))[1]})
rownames(TF_act) = TF_act$Cell
TF_act$Cell = NULL

load(paste0("post_traj_output_",cancer,"_all.RData"))

rows.in.both = intersect(rownames(TF_act),rownames(df))

#filter TF_act to include only cells in df
TF_act = TF_act[which(rownames(TF_act) %in% rows.in.both),]

#add pseudotime values to TF_act
TF_act$slingPseudotime_1 = sapply(rownames(TF_act), function(x) {df$slingPseudotime_1[which(rownames(df)==x)]})
TF_act$slingPseudotime_2 = sapply(rownames(TF_act), function(x) {df$slingPseudotime_2[which(rownames(df)==x)]})

#add DR values to TF_act
TF_act$DR1 = sapply(rownames(TF_act), function(x) {df$DR1[which(rownames(df)==x)]})
TF_act$DR2 = sapply(rownames(TF_act), function(x) {df$DR2[which(rownames(df)==x)]})

inds.TFs = 1:(ncol(TF_act)-4)


#compute TFs which correlate well with each trajectory
#iterate over 2 pseudotimes
for(i in 1:2){
  corr_ps = lapply(inds.TFs, function(x) {
    res = cor.test(TF_act[,x],TF_act[[paste0("slingPseudotime_",i)]], 
                   use = "complete.obs", method = "pearson")
    return(c(res$estimate, res$p.value))
  })
  df.corr = do.call("rbind",corr_ps)
  df.corr = as.data.frame(df.corr)
  rownames(df.corr) = colnames(TF_act)[inds.TFs]
  colnames(df.corr) = c("corr_ps","pval_1")
  
  Nplot = 5
  for(direc in c("up","down")){
    plot.list = list()
    plot.list2 = list()
    #plot top 5, both expression on trajectory and dot plot of expr v pseudotime
    for(b in 1:Nplot){
      if(direc=="up"){
        TF = rownames(df.corr[order(df.corr$corr_ps,decreasing=T),])[b]
      }else{
        TF = rownames(df.corr[order(df.corr$corr_ps),])[b]
      }
      plot.list2[[b]] = ggplot(TF_act, aes(x = .data[[paste0("slingPseudotime_",i)]], y = .data[[TF]])) +
        geom_point(size=0.001) + labs(title=TF)
      

      plot.list[[b]] = ggplot(TF_act, aes(x = DR1, y = DR2, color = .data[[TF]])) + geom_point(size=0.001)
    }
    p1 = wrap_plots(c(plot.list,plot.list2), ncol = Nplot, byrow = T)
    print(p1)
    ggsave(paste0(fig.dir,"/ps",i,"_TF_act_",direc,"_",cancer,".pdf"),
           width = 25, height = 5)
  }
}






