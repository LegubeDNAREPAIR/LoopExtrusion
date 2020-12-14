require(tidyverse)
require(plyranges)
require(rtracklayer)
require(BSgenome.Hsapiens.UCSC.hg19.masked)
seqlens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19.masked)
source("src/functions.R")
bless80 <- "/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed" %>% read_bed()

my_window <- 10000000
mes_seqlens <- seqlens %>% enframe()
mon_bed <- bless80 %>%
    anchor_center() %>% mutate(width = my_window*2) %>%
    as_tibble() %>% left_join(mes_seqlens,by = c("seqnames"="name")) %>% 
    filter(end <= value) %>% 
    as_granges() %>% anchor_center() %>% mutate(width=1)
mon_bed_custom <- mon_bed %>% sort()  %>% anchor_center() %>% mutate(width=2000) %>% mutate(DSB = name) %>% as_tibble() %>% dplyr::select(-value)%>% mutate(Group = "DSB") %>% as_granges()


damagedTAD_divided <- GRangesList(
    
    mon_bed %>% anchor_start() %>% stretch(9999) %>% mutate(Group = "DSB-10kb"),
    mon_bed %>% anchor_end() %>% stretch(9999) %>% mutate(Group = "DSB-10kb"),
    mon_bed %>% shift_right(10000) %>% anchor_start() %>% mutate(width=10000) %>% mutate(Group = "10kb-20kb"),
    mon_bed %>% shift_left(10000) %>% anchor_end() %>% mutate(width=10000) %>% mutate(Group = "10kb-20kb"),
    mon_bed %>% shift_right(30000) %>% anchor_start() %>% mutate(width=70000) %>% mutate(Group = "30kb-100kb"),
    mon_bed %>% shift_left(30000) %>% anchor_end() %>% mutate(width=70000) %>% mutate(Group = "30kb-100kb"),
    mon_bed %>% shift_right(100000) %>% anchor_start() %>% mutate(width=100000) %>% mutate(Group = "100kb-200kb"),
    mon_bed %>% shift_left(100000) %>% anchor_end() %>% mutate(width=100000) %>% mutate(Group = "100kb-200kb"),
    mon_bed %>% shift_right(200000) %>% anchor_start() %>% mutate(width=800000) %>% mutate(Group = "200kb-1Mb"),
    mon_bed %>% shift_left(200000) %>% anchor_end() %>% mutate(width=800000) %>% mutate(Group = "200kb-1Mb"),
    mon_bed %>% shift_right(1000000) %>% anchor_start() %>% mutate(width=4000000) %>% mutate(Group = "1Mb-5mb"),
    mon_bed %>% shift_left(1000000) %>% anchor_end() %>% mutate(width=4000000) %>% mutate(Group = "1Mb-5mb")
) %>% unlist() %>% sort() %>% as_tibble() %>% dplyr::select(seqnames,start,end,name,Group) %>% dplyr::rename(DSB = name) %>% as_granges()


#Load loop_anchor
mes_loops_files <- list.files("/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Hi-C/HiC/loop/hiccups/Damaged_Undamaged",full.names = T)
loop_data <- mes_loops_files %>% map(read_tsv,col_types = list(X1 = col_character(),X4 = col_character()),col_names = F)
names(loop_data) <- mes_loops_files %>% map(basename) %>% map(str_remove,"\\.bed")
mon_bed_custom <- mon_bed_custom %>% as_tibble() %>% select(-score) %>% as_granges()
loop_bed <- loop_data %>% map(process_loop)
loop_bed_inside <- loop_bed %>% map(function(x){do.call("c",x)}) %>% map(sort)%>% map(function(x){x %>% join_overlap_inner(damagedTAD_divided)})

#All loop size (not splided)
mes_bw <- c("/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Clouaire_HLNKYBGXC_SCC1//PROCESSED/mapping/EXPERIMENT/BIGWIG/HLNKYBGXC_Pool_ChIP-seq_legube_19s004478-1-1_Clouaire_lane1Rad21DIvA_sequence.exp_spikeinfactor.bw",
            "/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Clouaire_HLNKYBGXC_SCC1//PROCESSED/mapping/EXPERIMENT/BIGWIG/HLNKYBGXC_Pool_ChIP-seq_legube_19s004478-1-1_Clouaire_lane1Rad21OHT_sequence.exp_spikeinfactor.bw")
names(mes_bw) <- c("Rad21_DIvA","Rad21_OHT")
mes_bw <- mes_bw %>% map(import.bw,as="RleList")
res.loop <- lapply(loop_bed_inside,function(one_loop){
    one_loop <- one_loop %>% anchor_center() %>% mutate(width = 2000)
    lapply(names(mes_bw),function(one.w){
        message(one.w)
        rbind(
            Get1val(one.w,mes_bw[[one.w]],one_loop) %>% left_join(as_tibble(one_loop),by=c("rowname"="name")),
            Get1val(one.w,mes_bw[[one.w]],mon_bed_custom) %>% left_join(as_tibble(mon_bed_custom),by=c("rowname"="name")) %>% mutate(sizeGroup = ">2e+05")
        )
        
    }) %>% bind_rows()
}) %>% bind_rows(.id = "Type")
p3 <-res.loop %>%
    separate(wig,into = c("wig","Condition"))  %>%
    filter(str_detect(wig,"Rad21")) %>% 
    mutate(Group = fct_relevel(Group,c("DSB","DSB-10kb","10kb-20kb","30kb-100kb", "100kb-200kb","200kb-1Mb","1Mb-5mb"))) %>% 
    ggplot(aes(x=Group,y=value,fill = Condition)) + geom_boxplot() + theme_classic() +
    scale_y_log10()+
    facet_wrap(~wig,ncol=1,scales="free_y")
print(p3)
