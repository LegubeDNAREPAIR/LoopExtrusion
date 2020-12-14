.libPaths(c(.libPaths(),"/home/rochevin/R-packages/"))
setwd(snakemake@params[["dir"]]) ;
require(rtracklayer)
require(tidyverse)
require(BSgenome.Hsapiens.UCSC.hg19)
seqlens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
viewpoint <- snakemake@params[["viewpoint"]]
region.to.overlap <- GRanges(viewpoint%>% str_replace("_",":"))
inpath <- snakemake@params[["inpath"]]
extension <- snakemake@params[["extension"]] %>% as.numeric()
outpath <- snakemake@output[[1]]


bed.graph.files <- list.files(inpath,pattern = paste(viewpoint,"count",sep="."),full.names = T)
matrix.counts <- lapply(bed.graph.files,function(i){
	read_tsv(i,col_names = F) %>%
	rename(seqnames = X1) %>%
	rename(start = X2) %>%
	rename(end = X3) %>%
	rename(score = X7) %>% 
	mutate(pos = paste(seqnames,paste(start,end,sep="-"),sep=":")) %>%
	dplyr::select(pos,score)
})

#Merge matrix
matrix.counts.merge <- Reduce(function(x, y) {full_join(x, y,by="pos")}, matrix.counts) %>% tibble::column_to_rownames("pos")
colnames(matrix.counts.merge) <- bed.graph.files
matrix.counts.merge[is.na(matrix.counts.merge)] <- 0 #replace NA by 0 count
#EXPORT

write.table(matrix.counts.merge,file = outpath,quote=F,sep="\t",row.names = T,col.names = T)
