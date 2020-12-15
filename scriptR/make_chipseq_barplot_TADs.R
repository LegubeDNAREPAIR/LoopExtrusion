require(tidyverse)
require(rtracklayer)
require(plyranges)
source("src/functions.R")
TADb <- list.files("/media/HDD_ROCHER/PROJET_INGE/HiC_Coline",full.names = T,pattern="TADborder_.+\\.bed") %>% map(import.bed) %>% map(sort)
DSB174 <- read_bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/ASIsites_hg19_174_clived_IQ1.5.bed")



mes_bw <- PhDfunc::GetBWList()[c("MDC1_pOHT","GAM_clean_normalized_01022018_pOHT","53BP1_pOHT")]%>% str_replace("NAS","NAS1") %>% map(import.bw,as="RleList")
names(mes_bw) <- c("MDC1_pOHT","GAM_clean_normalized_01022018_pOHT","53BP1_pOHT")

TADs.DSB <- lapply(TADb,function(one_tad){
    
    tot.res <- GRanges()
    for(i in 1:length(DSB174)){
        res.p <- follow(DSB174[i],one_tad)
        if(is.na(res.p))
            next;
        
        pos.p <- tryCatch({
            one_tad[c(res.p,res.p-1,res.p-2,res.p-3,res.p-4)] %>% filter(seqnames == as.character(seqnames(DSB174[i]))) %>% start()
        },error = function(error){
            NULL;
        })
        pos.f <- tryCatch({
            one_tad[c(res.p+1,res.p+2,res.p+3,res.p+4,res.p+5)]  %>% filter(seqnames == as.character(seqnames(DSB174[i]))) %>% start()
        },error = function(error){
            NULL;
        })
        if(is.null(pos.p) | is.null(pos.f))
            next;
        
        if(length(pos.p)<5 | length(pos.f)<5)
            next;
        
        sres <- tibble(seqnames = as.character(seqnames(DSB174[i])),
                       start = c(
                           pos.p[5],
                           pos.p[4],
                           pos.p[3],
                           pos.p[2],
                           pos.p[1],
                           pos.f[1],
                           pos.f[2],
                           pos.f[3],
                           pos.f[4]
                       ),
                       end = c(
                           pos.p[4],
                           pos.p[3],
                           pos.p[2],
                           pos.p[1],
                           pos.f[1],
                           pos.f[2],
                           pos.f[3],
                           pos.f[4],
                           pos.f[5]
                       ),
                       name = str_c(DSB174[i]$name,c(
                           "TAD-4",
                           "TAD-3",
                           "TAD-2",
                           "TAD-1",
                           "TAD-Damaged",
                           "TAD+1",
                           "TAD+2",
                           "TAD+3",
                           "TAD+4"
                       ),sep="_")) %>% as_granges()
        tot.res <- c(tot.res,sres)
    }
    tot.res
}) %>% setNames(list.files("/media/HDD_ROCHER/PROJET_INGE/HiC_Coline",pattern="TADborder_.+\\.bed"))



cc_res <- lapply(names(mes_bw),function(one_n){
    one_w <- mes_bw[[one_n]]
    res <- lapply(names(TADs.DSB),function(one_b_n){
        one_bed <- TADs.DSB[[one_b_n]]
        Get1val(one_n,one_w,one_bed) %>% mutate(bed = one_b_n)
    }) %>% bind_rows()
}) %>% bind_rows()


#0.05 only
tadbn <- "TADborder_topdom_DIvA_manipA_hg19_50kb_ws10_strong_t-0.05.bed"
scc <- cc_res %>% filter(bed == tadbn) %>% mutate(TADpos = str_extract(rowname,"_TAD.+")) %>%
    mutate(TADpos = str_remove(TADpos,"_")) %>%
    mutate(rowname = str_remove(rowname,"_TAD.+")) %>% mutate(TADpos = factor(TADpos,levels = c("TAD-4","TAD-3","TAD-2","TAD-1","TAD-Damaged","TAD+1","TAD+2","TAD+3","TAD+4")))

# p1 <-scc %>% filter(wig == "53BP1_pOHT") %>% ggplot(aes(x=TADpos,y=value))+ geom_boxplot() + theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("53BP1_pOHT")
# p2 <-scc %>% filter(wig == "GAM_clean_normalized_01022018_pOHT") %>% ggplot(aes(x=TADpos,y=value))+ geom_boxplot() + theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("GAM_clean_normalized_01022018_pOHT")
# p3 <-scc %>% filter(wig == "MDC1_pOHT") %>% ggplot(aes(x=TADpos,y=value))+ geom_boxplot() + theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1))+ ggtitle("MDC1_pOHT")
p1 <-scc %>% filter(wig == "53BP1_pOHT") %>% ggplot(aes(x=TADpos,y=value)) + geom_bar(stat = "summary",fun.y="mean",col="black") + theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("53BP1_pOHT")
p2 <-scc %>% filter(wig == "GAM_clean_normalized_01022018_pOHT") %>% ggplot(aes(x=TADpos,y=value)) + geom_bar(stat = "summary",fun.y="mean",col="black") + theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("GAM_clean_normalized_01022018_pOHT")
p3 <-scc %>% filter(wig == "MDC1_pOHT") %>% ggplot(aes(x=TADpos,y=value)) + geom_bar(stat = "summary",fun.y="mean",col="black") + theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1))+ ggtitle("MDC1_pOHT")
