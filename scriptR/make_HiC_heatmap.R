require(HiTC)
require(tidyverse)
require(plyranges)
require(rtracklayer)
require(keras) # Only for %<-%
require(reshape2)
#Functions
source("src/functions.R")


obserOE = "OE" # observed, OE
if(obserOE == "OE"){
    manipA <- "rep1"
    manipB <- "rep2"
    my.lim <- 0.05
}else{
    manipA <- "manipA"
    manipB <- "manipB"
    my.lim <- 0.1
}

my_bed <- read_bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")

m.w <- "50kb"
sw <- 50000
my_lim <- 5000000
bin <- (my_lim/sw)+1

if(obserOE == "OE"){
    HiC.files <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/data/HTC/OE",pattern=m.w,full.names=T)
}else{
    HiC.files <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/data/HTC",pattern=m.w,full.names=T)
}

res.DSB <- lapply(unique(seqnames(my_bed)),function(chrom){
    message(chrom)
    f <- HiC.files[str_detect(HiC.files,str_c("_",chrom,"_"))]
    Replicate <- f %>% map(.%>% basename() %>% str_extract("manip[A-Z]+|rep[0-9]+"))
    diva.chr <- my_bed %>% filter(seqnames ==chrom) %>%
        anchor_center() %>% mutate(width = my_lim)
    
    if(length(diva.chr)<2)
        return(NULL)
    
    res <- lapply(unique(Replicate),function(j){
        s <- f[grep(j,f)]
        if(length(s)<2){
            NULL
        }else{
            message(s)
            Cond <- s %>% map(.%>% basename() %>% str_extract("DIvA|OHT"))
            c(res1,val1) %<-% get_sub_matrix_Hic(file = s[[1]],bed = diva.chr,binsize = bin)
            c(res2,val2) %<-% get_sub_matrix_Hic(file = s[[2]],bed = diva.chr,binsize = bin)
            
            # res1 <- lapply(res1,function(x){
            #     x * (val1/val2)
            # })
            
            res2 <- lapply(res2,function(x){
                x * (val1/val2)
            })
            
            res1 <- Reduce("+",res1)
            res2 <- Reduce("+",res2)
            full <- list(res1,res2)
            names(full) <- Cond
            full
        }
        
        
    })
    names(res) <- unique(Replicate)
    res
    
    
    
})


pOHT.A <- lapply(res.DSB,function(x){
    x[[manipA]][["OHT"]]
}) %>% plyr::compact()
pOHT.A <- Reduce("+",pOHT.A)
pOHT.B <- lapply(res.DSB,function(x){
    x[[manipB]][["OHT"]]
})%>% plyr::compact()
pOHT.B <- Reduce("+",pOHT.B)

mOHT.A <- lapply(res.DSB,function(x){
    x[[manipA]][["DIvA"]]
})%>% plyr::compact()
mOHT.A <- Reduce("+",mOHT.A)
mOHT.B <- lapply(res.DSB,function(x){
    x[[manipB]][["DIvA"]]
})%>% plyr::compact()
mOHT.B <- Reduce("+",mOHT.B)

data.plot <- list(
    "pOHT_A" = pOHT.A %>% as.matrix %>% melt,
    "mOHT_A" = mOHT.A %>% as.matrix %>% melt,
    "pOHT_B" = pOHT.B %>% as.matrix %>% melt,
    "mOHT_B" = mOHT.B %>% as.matrix %>% melt
) %>% bind_rows(.id = "Condition") %>%
    separate(Condition,into = c("Condition","Replicate"),sep = "_")%>%
    mutate(Var1 = forcats::lvls_revalue(Var1, as.character(1:bin))) %>%
    mutate(Var2 = forcats::lvls_revalue(Var2, as.character(1:bin)))

data.mean <- data.plot %>% spread(key = Replicate,value = value) %>%
    rowwise() %>%
    mutate(value = mean(A,B)) %>%
    dplyr::select(-A:-B)

data.ratio.AB <- data.plot %>% spread(key = Replicate,value = value) %>%
    rowwise() %>%
    mutate(meanValue = mean(c(A,B))) %>%
    dplyr::select(-A:-B) %>%
    spread(key = Condition,value = meanValue) %>%
    mutate(ratio = log2(pOHT/mOHT))

data.ratio <- data.plot %>%
    spread(key = Condition,value = value)%>%
    mutate(ratio = log2(pOHT/mOHT))
#Our plots comes from observed data
p.ratio.AB <- plot_ratio(data.ratio.AB,facet=FALSE,window = my_lim/2,my.quantile = 0.995,fixed = T,fixed.limite = c(my.lim,my.lim)) 
p.ratio <- plot_ratio(data.ratio,facet=TRUE,window = my_lim/2,my.quantile = 0.995,fixed = T,fixed.limite = c(my.lim,my.lim)) 
#Work only with OE data
p.n <- plot_telquel(data.plot,facet=TRUE,window = my_lim/2,my.quantile = 0.95,fixed = T,fixed.limite = 200) + scale_fill_gradientn(colors = RColorBrewer::brewer.pal(6, "YlOrRd"))
p.nmean <- plot_telquel(data.mean,facet=F,window = my_lim/2,my.quantile = 0.95,fixed = T,fixed.limite = 200) + scale_fill_gradientn(colors = RColorBrewer::brewer.pal(6, "YlOrRd")) + facet_wrap(~Condition)
