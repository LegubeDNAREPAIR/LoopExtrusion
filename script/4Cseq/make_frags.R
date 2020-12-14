.libPaths(c(.libPaths(),"/home/rochevin/R-packages/"))
setwd(snakemake@params[["dir"]]) ;
require(BSgenome.Hsapiens.UCSC.hg19)
require(Biostrings)
require(rtracklayer)
require(dplyr)




hg19 <- BSgenome.Hsapiens.UCSC.hg19
seqlens <- seqlengths( hg19 );
chromosomes <- seqnames(hg19)

MboI <- DNAString("GATC")
NlaIII <-DNAString("CATG")

match.mboI.NlaIII <- lapply(seqnames(hg19),function(chr){
	match.mboI <- matchPattern(MboI, hg19[[chr]], max.mismatch=0)
	match.NlaIII <- matchPattern(NlaIII, hg19[[chr]], max.mismatch=0)

	bed.mboI <- GRanges(seqnames = chr,ranges = match.mboI %>% ranges())
	bed.mboI$name <- "MboI"
	bed.NlaIII <- GRanges(seqnames = chr,ranges = match.NlaIII %>% ranges())
	bed.NlaIII$name <- "NlaIII"
	c(bed.mboI,bed.NlaIII) %>% sortSeqlevels() %>% sort()
}) %>% GRangesList() 


frags <- lapply(match.mboI.NlaIII,function(elm){
	frag <- GRanges(seqnames = seqnames(elm)@values,
				IRanges(
					start = end(elm[1:(length(elm)-1)]),
					end = start(elm[2:length(elm)])
				),
				name1 = elm[1:(length(elm)-1)]$name,
				name2 = elm[2:length(elm)]$name)
	#On exclut les fragments <80
	frag <- frag[which(width(frag)>80)]
	#On considère les MboI/NlaIII NlaIII/MboI
	frag$frag_type <- "bad"
	frag[which(frag$name1 != frag$name2)]$frag_type <- "good"

	#On exclut les positions des motifs
	start(frag) <- start(frag) +1
	end(frag) <- end(frag) -1
	frag
}) %>% do.call("c",.)

export.bed(frags[which(frags$frag_type == "good")],snakemake@output[[1]])
export.bed(frags[which(frags$frag_type == "bad")],snakemake@output[[2]])


