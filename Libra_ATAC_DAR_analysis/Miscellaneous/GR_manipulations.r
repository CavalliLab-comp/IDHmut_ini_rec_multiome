#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR4/bin/R
library(GenomicRanges)

# create a bedfile to the given path from a genomic
# ranges GRanges object
writeBedfromGR <- function(gr, path){
  df <- data.frame(seqnames = seqnames(gr),
    starts = start(gr)-1,
    ends = end(gr),
    names = c(rep(".", length(gr))),
    scores = c(rep(".", length(gr))),
    strands = strand(gr))

  write.table(df, file = path, quote = F, sep="\t", row.names = F, col.names=F)
}