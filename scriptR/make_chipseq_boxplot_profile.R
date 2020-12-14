require(dplyr)
require(rtracklayer)
require(ggplot2)
require(plyranges)
library("BSgenome.Hsapiens.UCSC.hg19")
seqlens = seqlengths( Hsapiens );
source("src/functions.R")
#Load AsiSI locations
asi = import.bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/ASIsites_hg19.bed")
bless80 = import.bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
HR = import.bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/BLESS_HR_JunFragPE_Rmdups_pm500bp.bed")
NHEJ = import.bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed")
Random80 = import.bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/80random.bed")
Random30 = import.bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/30random.bed")
uncut <- asi[!asi$name %in% bless80$name]

list.sites <- list("cut"=bless80,"uncut"=uncut,"HR"=HR,"NHEJ"=NHEJ,"Random80"=Random80,"Random30"=Random30)

#Get information from snakemake
window <- 2000
span <- 5

wigs <- lapply(c("/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Clouaire_HLNKYBGXC_SCC1//PROCESSED/mapping/EXPERIMENT/BIGWIG/HLNKYBGXC_Pool_ChIP-seq_legube_19s004478-1-1_Clouaire_lane1Rad21DIvA_sequence.exp_spikeinfactor.bw",
                 "/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Clouaire_HLNKYBGXC_SCC1//PROCESSED/mapping/EXPERIMENT/BIGWIG/HLNKYBGXC_Pool_ChIP-seq_legube_19s004478-1-1_Clouaire_lane1Rad21OHT_sequence.exp_spikeinfactor.bw"),import.bw,as="RleList")
names(wigs) <- c("Rad21_DIvA","Rad21_OHT")

dat.boxplot.all <- mclapply(wigs,function(wig){
    lapply(names(list.sites), ParaleliseViewboxplot,one.w=wig,list.sites=list.sites) %>%bind_rows()
},mc.cores = length(wigs)) %>% bind_rows(.id = "Condition")
prof.dat <- mclapply(wigs,function(wig){
    lapply(names(list.sites), ParaleliseViewprofile,one.w=wig,list.sites=list.sites) %>%bind_rows()
},mc.cores = length(wigs)) %>% bind_rows(.id = "Condition")





cutcolors <- c("#FDBECD","black","#BEBEBE")
HRNHEJ = c("#F5AB35","#049372","#BEBEBE")

filename <- "RAD21"
windowname <- "4kb"
#cutvsuncutvsrandom
p1 <- prof.dat %>% filter(Type %in% c("cut","uncut","Random80"))%>%ggplot(aes(Windows,Value,colour = Type)) +
    labs(list(title = "", x = "", y = "")) +
    geom_line()+
    facet_wrap(~ Condition,ncol = 1,scales = "free_x") +
    scale_colour_manual(values=cutcolors) +
    theme_classic() +
    ggtitle(paste(filename,windowname,sep="_"))
b1 <- dat.boxplot.all %>% filter(Type %in% c("cut","uncut","Random80"))%>%ggplot(aes(Type,score,fill = Type)) +
    labs(list(title = "", x = "", y = "")) +
    geom_boxplot()+
    facet_wrap(~ Condition,ncol = 2) +
    scale_fill_manual(values=cutcolors) +
    theme_classic() +
    ggtitle(paste(filename,windowname,sep="_"))
##HRvsNHEJ
p2 <- prof.dat %>% filter(Type %in% c("HR","NHEJ","Random30"))%>%ggplot(aes(Windows,Value,colour = Type)) +
    labs(list(title = "", x = "", y = "")) +
    geom_line()+
    facet_wrap(~ Condition,ncol = 1,scales = "free_x") +
    scale_colour_manual(values=HRNHEJ) +
    theme_classic() +
    ggtitle(paste(filename,windowname,sep="_"))
b2 <- dat.boxplot.all %>% filter(Type %in% c("HR","NHEJ","Random30"))%>%ggplot(aes(Type,score,fill = Type)) +
    labs(list(title = "", x = "", y = "")) +
    geom_boxplot()+
    facet_wrap(~ Condition,ncol = 2) +
    scale_fill_manual(values=HRNHEJ) +
    theme_classic() +
    ggtitle(paste(filename,windowname,sep="_"))

print(p1)
print(b1)
print(p2)
print(b2)

