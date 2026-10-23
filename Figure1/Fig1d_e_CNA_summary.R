# inferCNV for ATAC
#based on code from https://github.com/kokonech/mbOncoAberrations/blob/main/snATACseq/processAtacTumorSampleSV4.R

###########
## make figures of CNA called by inferCNV on ATAC data for each subtype
## use 200 malignant cells per sample
###########

rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
cancer = as.character(args[1])


library(Seurat)
library(Signac)
library(infercnv, lib.loc="/mnt/beegfs/home/alyne/persistent/R_libs")
library(data.table)
library(GenomicRanges)


Nss = 200 #number of malignant cells to take from each sample

#reda in list of all samples
samps.all = sort(read.table("~/persistent/oligo_astro_samples_with_malig.txt", header=F)[,2])

#read in list of samples which have ATAC
samps = fread(paste0("/mnt/beegfs/home/alyne/stagein/scenicplus_",cancer,"/sample_list.txt"),
              data.table=F, header=F)[,1]

setwd(paste0("/mnt/beegfs/home/alyne/stagein/data_tmp/numbat_inputs/infercnv_atac_",cancer))

#read in reference and malignant cells from each sample
raw.counts.list = list()
for(i in 1:length(samps)){
  raw.counts = ReadMtx(
    paste0(samps[i],"_peak_counts.mtx"),
    paste0(samps[i],"_peak_counts_cells.tsv"),
    paste0(samps[i],"_peak_counts_peaks.tsv"),
    cell.column = 1, feature.column = 1
  )
  raw.counts.list[[i]] = raw.counts[,sample(ncol(raw.counts),min(Nss,ncol(raw.counts)))]
  colnames(raw.counts.list[[i]]) = paste0(colnames(raw.counts.list[[i]]),"_",samps[i])
}
raw.counts = do.call("cbind",raw.counts.list)


ref.counts = ReadMtx(
  paste0("infercnv_atac_",cancer,"_reference_10000.mtx"),
  paste0("infercnv_atac_",cancer,"_reference_10000_cells.tsv"),
  paste0("infercnv_atac_",cancer,"_reference_10000_peaks.tsv"),
  cell.column = 1, feature.column = 1
)

#get raw counts for peaks in standard chromosomes
rawCounts <- cbind(raw.counts,ref.counts)
rawCounts <- rawCounts[ grep("chr[0-9]",rownames(rawCounts)),]


#make full annotations file by combining reference and sample-specific
annTable.ref = fread(paste0("infercnv_atac_",cancer,"_reference_10000_cell_annot.tsv"), data.table=F, header=F)
annTable.obs = data.frame(V1=colnames(raw.counts), V2 = paste0("malignant_",unlist(lapply(samps,function(x){rep(x,Nss)}))) )
annTable = rbind(annTable.ref,annTable.obs)
write.table(annTable, "all_samps_ATAC_cnv_ann.txt" ,col.names = F, row.names = F, sep="\t", quote=F)

#get reference cells as ones who don't have "malignant" st start of name
ref.cells = unique(annTable[,2])[!grepl("^malignant_",unique(annTable[,2]))]

#run infercnv
infercnv_obj = CreateInfercnvObject(rawCounts,
                                    annotations_file="all_samps_ATAC_cnv_ann.txt",
                                    delim="\t",
                                    gene_order_file=paste0("gene_order_",cancer,".tsv"),
                                    ref_group_names=ref.cells,
                                    chr_exclude=c("X","Y","MT") # default
)

out.dir = paste0("/mnt/beegfs/home/alyne/stageout/infercnv_outputs/infercnv_ATAC_reference_all_samps_",cancer)
dir.create(out.dir)

options(scipen = 100)
infercnv_obj = infercnv::run(infercnv_obj,
                             # cutoff: 1 for SmartSeq, 0.1 for 10x Genomics
                             cutoff=0.1,
                             out_dir=out.dir,
                             cluster_by_groups=TRUE,
                             cluster_references=FALSE,
                             k_obs_groups =  5,
                             smooth_method="runmeans",
                             analysis_mode="samples",
                             output_format = "pdf", # issue with atac
                             denoise=T,
                             plot_steps=F,
                             HMM=F,
                             write_expr_matrix = TRUE,
                             resume_mode=TRUE,
                             verbose=T
)


#replot figure with correct order for chromosomes
obj <- readRDS(paste0(out.dir,"/run.final.infercnv_obj"))
new_gene_order = data.frame()
for (chr_name in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) {
  new_gene_order = rbind(new_gene_order, obj@gene_order[which(obj@gene_order[["chr"]] == chr_name) , , drop=FALSE])
}
names(new_gene_order) <- c("chr", "start", "stop")
copy_infercnv_obj = obj
copy_infercnv_obj@gene_order = new_gene_order
copy_infercnv_obj@expr.data = obj@expr.data[rownames(new_gene_order), , drop=FALSE]


#plot in high resolution
plot_cnv(copy_infercnv_obj, out_dir=out.dir,
         cluster_by_groups=TRUE,
         cluster_references=FALSE,
         k_obs_groups =  5,
         contig_cex = 2.5,
         dynamic_resize = 0.01,
         output_filename = paste0("infercnv_all_samps_test_",cancer),
         output_format = "png",
         png_res = 1200,
         useRaster = TRUE)


