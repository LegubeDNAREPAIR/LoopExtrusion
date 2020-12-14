require(tidyverse)
require(cowplot)
require(plyranges)
require(rtracklayer)
require(BSgenome.Hsapiens.UCSC.hg19.masked)
require(ChIPseeker)
source("src/functions.R")

#RAD21/SMC1/SMC3 on loop anchor
bless80 <- "/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed" %>% read_bed()
my_window <- 1000000
mon_bed <- bless80 %>%
    anchor_center() %>% mutate(width = my_window*2)
mes_loops_files <- list.files("/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Hi-C/HiC/loop/hiccups/Damaged_Undamaged",full.names = T)
loop_data <- mes_loops_files %>% map(read_tsv,col_types = list(X1 = col_character(),X4 = col_character()),col_names = F)
names(loop_data) <- mes_loops_files %>% map(basename) %>% map(str_remove,"\\.bed")

loop_bed <- loop_data %>% map(process_loop)
mes_bw <- c("/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Clouaire_HLNKYBGXC_SCC1//PROCESSED/mapping/EXPERIMENT/BIGWIG/HLNKYBGXC_Pool_ChIP-seq_legube_19s004478-1-1_Clouaire_lane1Rad21DIvA_sequence.exp_spikeinfactor.bw",
            "/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Clouaire_HLNKYBGXC_SCC1//PROCESSED/mapping/EXPERIMENT/BIGWIG/HLNKYBGXC_Pool_ChIP-seq_legube_19s004478-1-1_Clouaire_lane1Rad21OHT_sequence.exp_spikeinfactor.bw",
            "/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/COHESIN/PhosphoCohesines/PROCESSED/mapping/BIGWIG/HGC5WBGXG_PhosphoMutant_Pool1_20s002497-1-1_Clouaire_lane1SMC1P966DIvA_sequence_normalized.bw",
            "/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/COHESIN/PhosphoCohesines/PROCESSED/mapping/BIGWIG/HGC5WBGXG_PhosphoMutant_Pool1_20s002497-1-1_Clouaire_lane1SMC1P966OHT_sequence_normalized.bw",
            "/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/COHESIN/PhosphoCohesines/PROCESSED/mapping/BIGWIG/HGC5WBGXG_PhosphoMutant_Pool1_20s002497-1-1_Clouaire_lane1SMC3P1083DIvA_sequence_normalized.bw",
            "/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/COHESIN/PhosphoCohesines/PROCESSED/mapping/BIGWIG/HGC5WBGXG_PhosphoMutant_Pool1_20s002497-1-1_Clouaire_lane1SMC3P1083OHT_sequence_normalized.bw")
names(mes_bw) <- basename(mes_bw) %>% 
    str_extract("SMC.+_sequence|Rad21DIvA|Rad21OHT") %>% str_remove("_sequence")
mes_bw <- mes_bw %>% map(import.bw,as="RleList")

loop_bed_verysmall <- loop_bed[c("loops_verysmall_DIvA_OHT_notdamaged_manipA","loops_verysmall_DIvA_OHT_damaged_manipA")] %>% map(function(x){do.call("c",x)}) %>% map(sort)
res.loop <- lapply(loop_bed_verysmall,function(one_loop){
    one_loop <- one_loop %>% anchor_center() %>% mutate(width = 2000)
    lapply(names(mes_bw),function(one.w){
        message(one.w)
        PhDfunc::Get1val(one.w,mes_bw[[one.w]],one_loop) %>% left_join(as_tibble(one_loop),by=c("rowname"="name"))
    }) %>% bind_rows()
}) %>% bind_rows(.id = "Type")
p.loop <- res.loop %>% 
    mutate(Condition = str_extract(wig,"DIvA|OHT")) %>% 
    mutate(wig = str_remove(wig,"DIvA|OHT")) %>% 
    spread(key = Condition,value = value) %>% mutate(ratio = log2(OHT/DIvA)) %>%
    ggplot(aes(x=Type,y=ratio)) + geom_boxplot() + theme_classic() + facet_wrap(~wig,ncol=1,scales="free_y") + geom_hline(yintercept = 0,col="red",linetype="dashed")
print(p.loop)



res.boxplot <- lapply(loop_bed_verysmall,function(one_loop){
    lapply(names(mes_bw),function(one.w){
        message(one.w)
        Get1val(one.w,mes_bw[[one.w]],one_loop) %>% left_join(as_tibble(one_loop),by=c("rowname"="name"))
    }) %>% bind_rows()
}) %>% bind_rows(.id = "Type")
res.boxplot <- res.boxplot %>% separate(wig,into=c("wig","Condition"),sep="_")

res.boxplot %>% 
    filter(wig == "Rad21") %>% 
    filter(sizeGroup == "<2e+05") %>% 
    mutate(Type = fct_relevel(Type,c("loops_verysmall_DIvA_OHT_damaged_manipA","loops_verysmall_DIvA_OHT_notdamaged_manipA"))) %>% 
    ggplot(aes(x=Condition,y=value)) +
    geom_boxplot(width=0.5,aes(fill=Condition)) + theme_classic() + facet_wrap(~Type,scales="free_x",nrow=1) + xlab("Type")
