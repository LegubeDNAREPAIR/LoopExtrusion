#Functions
#For HiC
loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}
get_sub_matrix_Hic <- function(file,bed,binsize){
    HTC <- loadRData(file)
    rows.chr <- x_intervals(HTC)
    val1 <- HTC %>% intdata %>% sum
    # HTC <- normPerReads(HTC)
    
    overlaps <- as(findOverlaps(bed,rows.chr), "List")
    overlaps.bin <- lapply(overlaps,length) %>% unlist()
    overlaps <- overlaps[overlaps.bin == binsize]
    
    res <- lapply(1:length(overlaps),function(i){
        x <- overlaps[[i]]
        if(as.character(strand(bed[i])) == "-"){
            x <- rev(x)
        }
        intdata(HTC)[x,x]
    })
    return(list(res,val1))
}
get_nb_read_by_chr <- function(file){
    HTC <- loadRData(file)
    HTC %>% intdata %>% sum
}
plot_ratio <-function(data.ratio,my.quantile = 0.95,facet=TRUE,window = 1000000,fixed = F,fixed.limite = c(0.1,0.1),DSB = "DSB"){
    if(fixed ==T){
        limite <- fixed.limite
        if(length(limite)<2){
            limite <- rep(limite,2)
        }
    }else{
        limite <- data.ratio %>% pull(ratio) %>% quantile(my.quantile) %>% as.numeric()
        limite <- rep(limite,2)
    }
    
    debut <- data.ratio$Var1 %>% levels %>% as.numeric() %>% min()
    milieu <- data.ratio$Var1 %>% levels %>% as.numeric() %>% median()
    fin <- data.ratio$Var1 %>% levels %>% as.numeric() %>% max()
    p.ratio <- data.ratio %>%
        mutate(ratio = ifelse(ratio > limite[2],limite[2],ratio)) %>%
        mutate(ratio = ifelse(ratio < -limite[1],-limite[1],ratio)) %>%
        ggplot(aes(x=Var1,y=Var2,fill=ratio)) + geom_tile() +
        scale_fill_gradient2(low = "#2980b9",high = "#f1c40f",midpoint = 0,mid = "black",limits = c(-limite[1],limite[2])) +
        scale_x_discrete(name = 'Position',
                         breaks = c(debut,milieu,fin),
                         labels = c(-window,DSB,window)
        ) +
        scale_y_discrete(name = 'Position',
                         breaks = c(debut,milieu,fin),
                         labels = c(-window,DSB,window)
        ) +
        # geom_vline(xintercept = milieu,linetype ="dashed",col="red")+
        theme_classic() +
        theme(legend.position="bottom")
    if(facet == TRUE){
        p.ratio + facet_wrap(~Replicate,ncol=1)
    }else{
        p.ratio
    }
}

plot_telquel <-function(data.plot,my.quantile = 0.90,facet=TRUE,window = 1000000,fixed = F,fixed.limite = 1000,DSB = "DSB"){
    if(fixed ==T){
        limite <- fixed.limite
        
    }else{
        limite <- data.plot %>% pull(value) %>% quantile(my.quantile) %>% as.numeric()
    }
    
    debut <- data.plot$Var1 %>% levels %>% as.numeric() %>% min()
    milieu <- data.plot$Var1 %>% levels %>% as.numeric() %>% median()
    fin <- data.plot$Var1 %>% levels %>% as.numeric() %>% max()
    p.ratio <- data.plot %>%
        mutate(value = ifelse(value > limite,limite,value)) %>%
        # mutate(value = log10(value)) %>%
        ggplot(aes(x=Var1,y=Var2,fill=value)) + geom_tile() +
        scale_fill_gradient2(low = "white",high = "#EA2027",limits = c(0,limite)) +
        scale_x_discrete(name = 'Position',
                         breaks = c(debut,milieu,fin),
                         labels = c(-window,DSB,window)
        ) +
        scale_y_discrete(name = 'Position',
                         breaks = c(debut,milieu,fin),
                         labels = c(-window,DSB,window)
        ) +
        # geom_vline(xintercept = milieu,linetype ="dashed",col="red")+
        theme_classic() +
        theme(legend.position="bottom")
    if(facet == TRUE){
        p.ratio + facet_wrap(~Replicate+Condition,ncol=2)
    }else{
        p.ratio
    }
}
#Boxplot 4CSeq
Get1val <- function (Name, one.w, x) 
{
    require(magrittr)
    lapply(split(x, droplevels(seqnames(x))), function(zz) {
        message(unique(as.character(seqnames(zz))))
        cov <- one.w[[unique(as.character(seqnames(zz)))]]
        score <- IRanges::Views(cov, start = start(zz), end = end(zz)) %>% 
            mean()
        tibble::tibble(wig = Name, value = score, rowname = zz$name)
    }) %>% bind_rows()
}
#Loop

process_loop <- function(my.table,sizeselector = 200000){
    x1 <- my.table %>% dplyr::select(X1,X2,X3) %>%
        dplyr::rename(seqnames = X1) %>%
        mutate(seqnames = str_c("chr",seqnames)) %>%
        dplyr::rename(start = X2) %>%
        dplyr::rename(end = X3) %>%
        as_granges() %>% mutate(name = str_c("chr1_",1:nrow(my.table)))
    x2 <- my.table %>% dplyr::select(X4,X5,X6) %>%
        dplyr::rename(seqnames = X4) %>%
        mutate(seqnames = str_c("chr",seqnames)) %>%
        dplyr::rename(start = X5) %>%
        dplyr::rename(end = X6) %>%
        as_granges() %>% mutate(name = str_c("chr2_",1:nrow(my.table)))
    sizeLoop <- abs(start(x1) - end(x2))
    sizeGroup <- ifelse(sizeLoop < sizeselector,str_c("<",sizeselector),str_c(">",sizeselector))
    x1$sizeGroup <- sizeGroup
    x2$sizeGroup <- sizeGroup
    return(list(x1,x2))
}
#BIGWIG
#process BW
computeProfile = function( bed, wig, w = 20000, span = 200, seqlens ,method="mean"){
    if( class( wig ) != "SimpleRleList" ){
        stop( "ERROR : unknown class of wig, please provide a SimpleRleList" );
    }
    mat = NULL;
    for( i in 1:length( bed ) ){
        message( i, "/", length( bed ) );
        bedi = bed[i, ];
        chr = as.character( seqnames( bedi ) );
        cov = wig[[chr]];
        
        
        center = start( bedi ) + 4;
        stW = center - w;
        edW = center + w;
        
        if( span == 1 ){
            vm = as.numeric( Views( cov, start = stW, end = edW )[[1]] )
        }else{
            sts = seq( stW, edW - span + 1, span );
            eds = seq( stW + span - 1, edW, span );
            v = Views( cov, start = sts, end = eds );
            if(method =="sum"){
                vm =  sum( v );    
            }else {
                vm =  mean( v );
            }
            
            vm[sts >= seqlens[chr] | sts > length( cov )] = 0;
            vm[sts < 1] = 0;
        }
        mat = rbind( mat, vm );
    }
    #rv = colMeans( mat );
    return( mat );
}


ParaleliseViewboxplot <- function(x,one.w,list.sites){
    subdom <- list.sites[[x]]
    lapply(split(subdom,droplevels(seqnames(subdom))),function(zz){
        zz <- zz %>% anchor_center() %>% mutate(width = 2*window)
        message(unique(as.character(seqnames(zz))))
        cov <- one.w[[unique(as.character(seqnames(zz)))]]
        scoredat <- Views( cov, start = start(zz), end = end(zz) )%>% sum()
        zz %>% as_tibble() %>% mutate(score = scoredat) %>% mutate(Type = x)
    }) %>% bind_rows()
}

ParaleliseViewprofile <- function(x,one.w,list.sites){
    subdom <- list.sites[[x]]
    d1 <- computeProfile(bed=subdom,wig = one.w,seqlens = seqlens,w = window,span=span) %>% colMeans()
    data.frame(Value=d1,Windows=seq(-window,window-span+1,span),Type=x)
}