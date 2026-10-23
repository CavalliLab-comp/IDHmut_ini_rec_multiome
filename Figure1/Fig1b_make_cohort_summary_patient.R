library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(readxl)
library(tidyr)
library(ggnewscale)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

setwd("/mnt/beegfs/home/alyne/stagein")
fig.dir = "/mnt/beegfs/home/alyne/stageout/figures_cohort_summary"

#read in sample info
samp.md <- as.data.frame(read_excel("../persistent/sample_summary_astro_oligo_anon.xlsx", col_types = "text"))
samp.md$DOB <- as.Date(as.numeric(samp.md$DOB), origin = "1899-12-30")
samp.md$Surgery_date <- as.Date(as.numeric(samp.md$Surgery_date), origin = "1899-12-30")
samp.md$DDN <- as.Date(as.numeric(samp.md$DDN), origin = "1899-12-30")


df = data.frame(patient = unique(samp.md$patient)[!is.na(unique(samp.md$patient))])


samps.rna.only = c("oligo_P0011_ini_S0022", "astro_P0017_ini_S0032", "astro_P0018_rec1_S0035", "astro_P0025_ini_S0048")
samp.no.malig = c("oligo_P0003_rec2_S0007","astro_P0014_ini_S0026","oligo_P0022_rec1_S0043","oligo_P0030_ini_S0058")
samp.not.work = c("oligo_P0003_ini_S0008","astro_P0026_rec1_S0051","astro_P0013_ini_S0025")


#make wording a bit cleaner for paper
df$sex = sapply(df$patient, function(x) {samp.md$sex[which(samp.md$patient==x)[1]]})
df$age = sapply(df$patient, function(x) {as.numeric(samp.md$Age[which(samp.md$patient==x)[1]])})
df$cancer = sapply(df$patient, function(x) {samp.md$cancer[which(samp.md$patient==x)[1]]})
df$cancer[which(df$cancer=="oligo")] = "Oligodendroglioma"
df$cancer[which(df$cancer=="astro")] = "Astrocytoma"
df$grade1 = sapply(df$patient, function(x) {
  if(length(which(samp.md$patient==x & samp.md$ini_rec=="ini"))==1){
    samp.md$Grade[which(samp.md$patient==x & samp.md$ini_rec=="ini")]
  }else{
    NA
  }
})
df$grade2 = sapply(df$patient, function(x) {
  if(length(which(samp.md$patient==x & samp.md$ini_rec=="rec1"))==1){
    samp.md$Grade[which(samp.md$patient==x & samp.md$ini_rec=="rec1")]
  }else{
    NA
  }
})
df$grade3 = sapply(df$patient, function(x) {
  if(length(which(samp.md$patient==x & samp.md$ini_rec=="rec2"))==1){
    samp.md$Grade[which(samp.md$patient==x & samp.md$ini_rec=="rec2")]
  }else{
    NA
  }
})
df$multiome1 = sapply(df$patient, function(x) {
  if(length(which(samp.md$patient==x & samp.md$ini_rec=="ini"))==1){
    ind = which(samp.md$patient==x & samp.md$ini_rec=="ini")
    if(samp.md$Sample[ind] %in% samp.no.malig){
      "No malignant cells"
    }else{
      if(samp.md$Sample[ind] %in% samps.rna.only){
        "RNA only"
      }else{"Multiome"}
    }
  }else{NA}
})
df$multiome2 = sapply(df$patient, function(x) {
  if(length(which(samp.md$patient==x & samp.md$ini_rec=="rec1"))==1){
    ind = which(samp.md$patient==x & samp.md$ini_rec=="rec1")
    if(samp.md$Sample[ind] %in% samp.no.malig){
      "No malignant cells"
    }else{
      if(samp.md$Sample[ind] %in% samps.rna.only){
        "RNA only"
      }else{"Multiome"}
    }
  }else{NA}
})
df$multiome3 = sapply(df$patient, function(x) {
  if(length(which(samp.md$patient==x & samp.md$ini_rec=="rec2"))==1){
    ind = which(samp.md$patient==x & samp.md$ini_rec=="rec2")
    if(samp.md$Sample[ind] %in% samp.no.malig){
      "No malignant cells"
    }else{
      if(samp.md$Sample[ind] %in% samps.rna.only){
        "RNA only"
      }else{"Multiome"}
    }
  }else{NA}
})


#compute recurrence and final time points in years
df$time_ini = 0
df$time_rec1 = sapply(df$patient, function(x) {
  if(length(which(samp.md$patient==x & samp.md$ini_rec=="rec1"))==1){
    as.numeric(samp.md$Surgery_date[which(samp.md$patient==x & samp.md$ini_rec=="rec1")] - 
      samp.md$Surgery_date[which(samp.md$patient==x & samp.md$ini_rec=="ini")]) /365
  }else{
    NA
  }
})
df$time_rec2 = sapply(df$patient, function(x) {
  if(length(which(samp.md$patient==x & samp.md$ini_rec=="rec2"))==1){
    as.numeric(samp.md$Surgery_date[which(samp.md$patient==x & samp.md$ini_rec=="rec2")] - 
                 samp.md$Surgery_date[which(samp.md$patient==x & samp.md$ini_rec=="ini")]) /365
  }else{
    NA
  }
})
df$time_ddn = sapply(df$patient, function(x) {
  rel.rows = samp.md$DDN[which(samp.md$patient==x)]
  as.numeric(rel.rows[which(!is.na(rel.rows))] - samp.md$Surgery_date[which(samp.md$patient==x & samp.md$ini_rec=="ini")])/365
})
df$state_ddn = sapply(df$patient, function(x) {
  rel.rows = samp.md$State_DDN[which(samp.md$patient==x)]
  if(rel.rows[which(!is.na(rel.rows))] %in% c("Stable","Progression")){
    "Alive"
  }else if(rel.rows[which(!is.na(rel.rows))]=="Dead"){
    "Dead"
  }
})


#get adjuvant treatment status in between t0, t1, and t2
df$PCV1 = sapply(df$patient, function(x) {
  if(length(which(samp.md$patient==x & samp.md$ini_rec=="rec1"))==1){
    samp.md$PCV[which(samp.md$patient==x & samp.md$ini_rec=="rec1")]
  }else{NA}
})
df$TMZ1 = sapply(df$patient, function(x) {
  if(length(which(samp.md$patient==x & samp.md$ini_rec=="rec1"))==1){
    samp.md$TMZ[which(samp.md$patient==x & samp.md$ini_rec=="rec1")]
  }else{NA}
})
df$RT1 = sapply(df$patient, function(x) {
  if(length(which(samp.md$patient==x & samp.md$ini_rec=="rec1"))==1){
    samp.md$Radiotherapy[which(samp.md$patient==x & samp.md$ini_rec=="rec1")]
  }else{NA}
})
df$PCV2 = sapply(df$patient, function(x) {
  if(length(which(samp.md$patient==x & samp.md$ini_rec=="rec2"))==1){
    samp.md$PCV[which(samp.md$patient==x & samp.md$ini_rec=="rec2")]
  }else{NA}
})
df$TMZ2 = sapply(df$patient, function(x) {
  if(length(which(samp.md$patient==x & samp.md$ini_rec=="rec2"))==1){
    samp.md$TMZ[which(samp.md$patient==x & samp.md$ini_rec=="rec2")]
  }else{NA}
})
df$RT2 = sapply(df$patient, function(x) {
  if(length(which(samp.md$patient==x & samp.md$ini_rec=="rec2"))==1){
    samp.md$Radiotherapy[which(samp.md$patient==x & samp.md$ini_rec=="rec2")]
  }else{NA}
})
#remove samples which didn't work at all
samp.not.work = c("oligo_P0003_ini_S0008","astro_P0026_rec1_S0051","astro_P0013_ini_S0025")
df[df$patient=="P0003",c("multiome1")] = NA
df[df$patient=="P0026",c("multiome2")] = NA
df = df[df$patient!="P0013",]


#tile plot with cancer, patient, sex

#make colour scales
cols.mut = c("1" = "grey30", "0" = "gray95")
cols.grade = c("2" = "skyblue", "3" = "royalblue", "4" = "darkblue")
cols.pat = gg_color_hue(length(unique(df$patient)))
#order pats by cancer and ini_rec
pat.ord = unique(df$patient[order(df$cancer, df$patient)])
names(cols.pat) = c(matrix(pat.ord,ncol=2,byrow=T))
cols.ddnstate = c("Alive" = "coral", "Dead" = "aquamarine3")
cols.ini_rec = c("ini" = "violet", "rec1" = "purple", "rec2" = "purple4")
cols.cancer = c("Astrocytoma" = "darkgreen", "Oligodendroglioma" = "lightgreen")
cols.sex = c("F" = "red3", "M" = "blue3")
cols.data = c("Multiome" = "darkorange4", "RNA only" = "chocolate2", "No malignant cells" = "burlywood1")
full.color.scale = c(cols.mut, cols.grade, cols.pat, cols.ini_rec, cols.cancer, cols.sex, cols.ddnstate, cols.data)



####make plot


#first astrocytomas
df.oligo = df[df$cancer=="Astrocytoma",]
plot.df <- df.oligo |>
  mutate(patient_id = patient) |>
  pivot_longer(cols = c(patient, sex, cancer),names_to = "variable", values_to = "value")
plot.df$variable[which(plot.df$variable=="cancer")] = "subtype"
p1 = ggplot(plot.df, aes(x = patient_id, y = variable, fill = value)) +
  geom_tile(color = "white") + theme_minimal() + scale_fill_manual(values = full.color.scale) +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill="none")

#dot plot with age at first surgery
p2 = ggplot(df.oligo, aes(x = patient, y = age)) + geom_point(size = 3) + lims(y = c(25,75)) + labs(y="Age at initial surgery (yrs)") +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

#tile plot with multiome and treatment info
plot.treat <- pivot_longer(df.oligo, cols = c(multiome1, PCV1, TMZ1, RT1, multiome2, PCV2, TMZ2, RT2, multiome3),
                           names_to = "variable",values_to = "value")
plot.treat$variable = factor(plot.treat$variable, levels = c("multiome1", "PCV1", "TMZ1", "RT1", 
                                                             "multiome2", "PCV2", "TMZ2", "RT2", "multiome3"))
plot.treat$value = factor(plot.treat$value, levels = c("0","1","Multiome","RNA only","No malignant cells"))
p3 = ggplot(plot.treat, aes(x = patient, y = variable, fill = value)) + guides(fill="none") +
  geom_tile(color = "white") + theme_minimal() + scale_fill_manual(values = full.color.scale, na.value = "white") + labs(y="Treatment") +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + scale_y_discrete(
    labels = c( multiome1 = "Initial",PCV1 = "PCV",TMZ1 = "TMZ",RT1 = "RT",multiome2 = "Recurrence 1",PCV2 = "PCV",TMZ2 = "TMZ",RT2 = "RT",multiome3 = "Recurrence 2")
  )


#dot plot with recurrence times and DDN info, coloured by grade
plot.df2 <- bind_rows(
  transmute(df.oligo, patient, time = time_ini, event = "surgery", state = grade1),
  transmute(df.oligo, patient, time = time_rec1, event = "surgery", state = grade2),
  transmute(df.oligo, patient, time = time_rec2, event = "surgery", state = grade3),
  transmute(df.oligo, patient, time = time_ddn, event = "DDN", state = state_ddn)
)
plot.df2$event = factor(plot.df2$event, levels = c("surgery","DDN"))
p4 = ggplot(plot.df2, aes(x = patient, y = time, colour = state, shape = event)) + geom_point(size = 3) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + lims(y = c(0,28.1)) + labs(y="Time (yrs)") +
  scale_color_manual(values = full.color.scale) + scale_shape_manual(values = c(surgery = 16, DDN = 15))  + guides(color="none",shape="none")


#then oligodendrogliomas
df.oligo = df[df$cancer=="Oligodendroglioma",]
plot.df <- df.oligo |> mutate(patient_id = patient) |>
  pivot_longer(cols = c(patient, sex, cancer),names_to = "variable", values_to = "value")
plot.df$variable[which(plot.df$variable=="cancer")] = "subtype"
p5 = ggplot(plot.df, aes(x = patient_id, y = variable, fill = value)) +
  geom_tile(color = "white") + theme_minimal() + scale_fill_manual(values = full.color.scale) +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1),
        legend.text = element_text(size=6),legend.key.size = unit(0.3, "cm"),axis.title.y=element_blank(),
        axis.text.y=element_blank(),axis.ticks.y=element_blank())

#dot plot with age at first surgery
p6 = ggplot(df.oligo, aes(x = patient, y = age)) + geom_point(size = 3) + lims(y = c(25,75)) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        legend.text = element_text(size=6),legend.key.size = unit(0.3, "cm"),axis.title.y=element_blank(),
        axis.text.y=element_blank(),axis.ticks.y=element_blank())

#tile plot with multiome and treatment info
plot.treat <- pivot_longer(df.oligo, cols = c(multiome1, PCV1, TMZ1, RT1, multiome2, PCV2, TMZ2, RT2, multiome3),
                           names_to = "variable",values_to = "value")
plot.treat$variable = factor(plot.treat$variable, levels = c("multiome1", "PCV1", "TMZ1", "RT1", 
                                                             "multiome2", "PCV2", "TMZ2", "RT2", "multiome3"))
plot.treat$value = factor(plot.treat$value, levels = c("0","1","Multiome","RNA only","No malignant cells"))
p7 = ggplot(plot.treat, aes(x = patient, y = variable, fill = value)) +
  geom_tile(color = "white") + theme_minimal() + scale_fill_manual(values = full.color.scale, na.value = "white") +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        legend.text = element_text(size=6),legend.key.size = unit(0.3, "cm"),axis.title.y=element_blank(),
        axis.text.y=element_blank(),axis.ticks.y=element_blank())

#dot plot with recurrence times and DDN info, coloured by grade
plot.df2 <- bind_rows(
  transmute(df.oligo, patient, time = time_ini, event = "surgery", state = grade1),
  transmute(df.oligo, patient, time = time_rec1, event = "surgery", state = grade2),
  transmute(df.oligo, patient, time = time_rec2, event = "surgery", state = grade3),
  transmute(df.oligo, patient, time = time_ddn, event = "DDN", state = state_ddn)
)
plot.df2$event = factor(plot.df2$event, levels = c("surgery","DDN"))
p8 = ggplot(plot.df2, aes(x = patient, y = time, colour = state, shape = event)) + geom_point(size = 3) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        legend.text = element_text(size=6),legend.key.size = unit(0.3, "cm"),axis.title.y=element_blank(),
        axis.text.y=element_blank(),axis.ticks.y=element_blank()) + lims(y = c(0,28.1)) +
  scale_color_manual(values = full.color.scale) + scale_shape_manual(values = c(surgery = 16, DDN = 15))

wrap_plots(list(p4,p3,p2,p1,p8,p7,p6,p5), ncol=2, byrow = F, heights = c(5,3,3,1), widths = c(5,8))
ggsave(paste0(fig.dir,"/cohort_patient_summary.pdf"), width = 11)





