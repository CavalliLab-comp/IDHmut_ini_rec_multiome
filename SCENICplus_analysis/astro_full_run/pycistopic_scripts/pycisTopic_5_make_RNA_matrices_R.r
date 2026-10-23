library(Seurat)
library(Signac)
library(data.table)
library(future)
options(future.globals.maxSize = 100 * 1024 ^ 3)

cancer = "astro"

base.dir = paste0("/mnt/beegfs/home/alyne")
setwd(base.dir)


#read in samples
samples.newnames = read.table(paste0("stagein/scenicplus_",cancer,"/sample_list.txt"), header = F, 
                              colClasses = c("character"), stringsAsFactors = F)[,1]

#read in merged object
oligo.ss = readRDS(paste0("stagein/scenicplus_",cancer,"/oligo.ss_integrated_RPCA_all_",cancer,"_doub_rem.rds"))
# oligo.ss = JoinLayers(oligo.ss)


#find variable features
DefaultAssay(oligo.ss) = "RNA"

#number of genes expressed in 1% of cells
ncell = ncol(oligo.ss[['RNA']]$counts)
genes = rownames(oligo.ss[['RNA']]$counts)
mat = oligo.ss[['RNA']]$counts
genes_detect = 100*Matrix::rowSums(mat > 0)/ncell
genes_keep = rownames(mat)[which(genes_detect>1)]


#read in TFs
data.tfs.all = read.table("/mnt/beegfs/home/alyne/stagein/required_files/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl",
                          sep="\t", comment.char = "", header=T)

#find percentage of cells each TF expressed in
#find TFs in count matrix
tfs <- intersect(unique(data.tfs.all$gene_name), rownames(mat))
#compute percentage
tf_detect <- 100*Matrix::rowSums(mat[tfs, , drop = FALSE] > 0)/ncell


#keep all TFs, since aren't that many and I guess low expressed will be filtered anyway
feats.keep = union(genes_keep,tfs)


#extract matrices
to_write_sct = as.data.frame(t(as.matrix(oligo.ss[['SCT']]$data[rownames(oligo.ss[['SCT']]$data) %in% feats.keep,])))
to_write_rna = as.data.frame(t(as.matrix(oligo.ss[['RNA']]$counts[rownames(oligo.ss[['RNA']]$counts) %in% feats.keep,])))


#remove samples which only have RNA
to_write_sct = to_write_sct[which(oligo.ss$orig.ident %in% samples.newnames),]
to_write_rna = to_write_rna[which(oligo.ss$orig.ident %in% samples.newnames),]


#read in list of cells which passed ATAC QC
cells.ATAC = fread(paste0("stageout/pycisTopic_output_",cancer,"/cells_passing_qc_ATAC.csv"), data.table=F, 
                   header=T, na.strings="NaN")
cells.ATAC.rename = lapply(seq(ncol(cells.ATAC)), function(x) {
  paste(cells.ATAC[which(!is.na(cells.ATAC[,x])),x],colnames(cells.ATAC)[x],sep="-")
})


#filter cells in oligo.ss to keep only those passing pycistopic qc
to_write_sct = to_write_sct[which(rownames(to_write_sct) %in% unlist(cells.ATAC.rename)),]
to_write_rna = to_write_rna[which(rownames(to_write_rna) %in% unlist(cells.ATAC.rename)),]


#write array of final cells passing both my and pycistopic's QC
cell.list.QC = lapply(samples.newnames, function(x) {
  cell.names = names(oligo.ss$orig.ident)[which(oligo.ss$orig.ident==x)]
  cell.names = cell.names[which(cell.names %in% rownames(to_write_rna))]
  sapply(cell.names, function(y) {
    paste(unlist(strsplit(y,"-"))[1:2],collapse="-")
  })
})
names(cell.list.QC) = samples.newnames
n <- max(sapply(cell.list.QC, length))
mat.to.print = do.call(cbind, lapply(cell.list.QC, function(x) x[1:n, drop = FALSE]))


dir.RNA = paste0("stagein/scenicplus_",cancer,"/RNA")
dir.create(dir.RNA)
fwrite(mat.to.print,
       file = paste0(dir.RNA,"/cells_passing_both_QC_",cancer,".tsv.gz"),
       sep = "\t", append = F, col.names = T, row.names = F, quote = F)



#change names to match pycistopic
rownames(to_write_sct) = paste(rownames(to_write_sct), sapply(rownames(to_write_sct), function(x) {unlist(strsplit(x,"-"))[3]}),
                               sep="___")
fwrite(to_write_sct,
       file = paste0(dir.RNA,"/SCT_data_scenicplus_",cancer,".tsv.gz"),
       sep = "\t", append = F, col.names = T, row.names = T, quote = F)


rownames(to_write_rna) = paste(rownames(to_write_rna), sapply(rownames(to_write_rna), function(x) {unlist(strsplit(x,"-"))[3]}),
                               sep="___")
fwrite(to_write_rna,
       file = paste0(dir.RNA,"/raw_counts_scenicplus_",cancer,".tsv.gz"),
       sep = "\t", append = F, col.names = T, row.names = T, quote = F)












############
### old way using variable features
############

# oligo.ss = FindVariableFeatures(oligo.ss, nfeatures = 3500)
# var.feats = VariableFeatures(oligo.ss)
# 
# 
# #read in TFs
# data.tfs.all = read.table("/mnt/beegfs/home/alyne/stagein/required_files/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl",
#                           sep="\t", comment.char = "", header=T)
# 
# #find percentage of cells each TF expressed in
# ncell = ncol(oligo.ss[['SCT']]$data)
# mat = oligo.ss[['RNA']]$counts
# #find TFs in count matrix
# tfs <- intersect(unique(data.tfs.all$gene_name), rownames(mat))
# #compute percentage
# tf_detect <- 100*Matrix::rowSums(mat[tfs, , drop = FALSE] > 0)/ncell
# 
# 
# #add the ones detected in more than 2% cells
# var.feats = union(var.feats,names(tf_detect)[which(tf_detect>1)])





