smoothy <- function(x,values){
    vect <- values[["vect"]];
    rw <- values[["rw"]]
    res <- NULL;
    if(x<=rw){
        res <- mean(vect[1:(rw+x)]);
    }else if (x>(length(vect)-rw)){
        res <- mean(vect[(x-rw):length(vect)]);
    }else {
        res <- mean(vect[(x-rw):(x+rw)]);
    }
    return(res);
}

setwd(snakemake@params[["dir"]]) ;
library(Rsamtools) ; # ScanBam
library(GenomicAlignments) ; # GAlignement
library("BSgenome.Hsapiens.UCSC.hg19") ;
library( "rtracklayer" ) ; #Export bw
seqlens <- seqlengths( Hsapiens );
myWig <- import.bw(snakemake@input[[1]], as = "RleList");

w <- as.numeric(snakemake@params[["smooth_windows"]]) ;
sw <- as.numeric(snakemake@params[["smooth_subwindows"]]) ;
message(snakemake@input[[1]])
message(w)
message(sw)


wig <- GRanges() ;

for( chr in names(myWig) ){
    message( chr );
    cov = myWig[[chr]] ;
    len = length( cov );
    
    st <- seq( 1, len - sw + 1, sw ) ;
    ed <- seq( sw, len, sw ) ;
    v <- Views( cov, start = st, end = ed ) ;
    sc <- mean( v, na.rm = TRUE ) ;
    sc[is.nan( sc )] <- 0 ;
    
    res = sapply(1:length(sc),smoothy,values=list("vect"=sc,"rw"=(w/sw)))
    
    
    wig <- c( wig, GRanges( seqnames = chr, strand = "*", ranges = ranges( v ), score = res ) ) ;
}
seqlengths( wig ) <- seqlens[names(seqlengths( wig ) )] ;
export.bw( wig[!is.na(wig$score)], snakemake@output[[1]] )


