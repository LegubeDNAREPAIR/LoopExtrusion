require(BSgenome.Hsapiens.UCSC.hg19.masked)
require(regioneR)
require(tidyverse)
require(plyranges)


RAD21.narrow <- "/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Clouaire_HLNKYBGXC_SCC1/PROCESSED/mapping/EXPERIMENT/MACS/" %>%
    list.files(pattern=".narrowPeak",full.names = T) 
names(RAD21.narrow) <- RAD21.narrow %>% basename() %>% str_remove(".narrowPeak")
RAD21.narrow <- RAD21.narrow  %>% map(read_tsv,col_names = F) %>% bind_rows(.id = "TypePeak")
names(RAD21.narrow) <- c(
    "TypePeak",
    "seqnames",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "FC_PeakSummit",
    "-log10pvalue_PeakSummit",
    " -log10qvalue_PeakSummit",
    "RelativeDistanceToSummit"
)
RAD21.narrow$strand <- "*"


RAD21.narrow.gr <- RAD21.narrow %>% as_granges() 


RAD21.narrow.gr <- RAD21.narrow.gr %>% sort %>% filter(!str_detect(name,"Rad21OHT")) 


mes_bw <- lapply(c("/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Clouaire_HLNKYBGXC_SCC1//PROCESSED/mapping/EXPERIMENT/BIGWIG/HLNKYBGXC_Pool_ChIP-seq_legube_19s004478-1-1_Clouaire_lane1Rad21DIvA_sequence.exp_spikeinfactor.bw",
                 "/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Clouaire_HLNKYBGXC_SCC1//PROCESSED/mapping/EXPERIMENT/BIGWIG/HLNKYBGXC_Pool_ChIP-seq_legube_19s004478-1-1_Clouaire_lane1Rad21OHT_sequence.exp_spikeinfactor.bw"),import.bw,as="RleList")
names(mes_bw) <- c("Rad21_DIvA","Rad21_OHT")

res.boxplot.bg <- lapply(names(mes_bw),function(one.w){
    message(one.w)
    Get1val(one.w,mes_bw[[one.w]],RAD21.narrow.gr)
}) %>% bind_rows()


res.boxplot.bg <- res.boxplot.bg %>% separate(wig,into=c("wig","Condition"),sep="_") %>% 
    mutate(Peak = str_extract(rowname,"Rad21DIvA|interPeaks|backgroundRegions"))



#FILTER BY DAMAGED TADS
#WITHOUT RAD21pOHT peaks
damaged.tad <- read_bed("/mnt/NAS1/DATA/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed") %>% anchor_center() %>% mutate(width=1000000)

RAD21.mOHT.damaged.tad <- RAD21.narrow.gr %>% filter_by_overlaps(damaged.tad)
RAD21.mOHT.undamaged.tad <- RAD21.narrow.gr%>% filter_by_non_overlaps(damaged.tad) 

ccfilter <- .  %>% as_tibble() %>% dplyr::select(name)

res.boxplot.bg.tad <- list(
    "Peaks_damaged"=RAD21.mOHT.damaged.tad %>% ccfilter,
    "Peaks_undamaged"=RAD21.mOHT.undamaged.tad %>% ccfilter,
    "Peaks_all"=RAD21.narrow.gr %>% ccfilter
) %>% bind_rows(.id="Type") %>% inner_join(res.boxplot.bg,by=c("name"="rowname"))

res.boxplot.bg.tad <- res.boxplot.bg.tad %>% separate(Type,into = c("Peak","TAD")) 

p.bg.rad21.RAD21 <- res.boxplot.bg.tad%>% filter(Peak =="Peaks") %>% 
    ggplot(aes(x=Condition,y=value)) +
    theme_classic()+facet_grid(TAD~Peak) + theme(legend.position="none")

print(p.bg.rad21.RAD21 + geom_violin(aes(col=Condition)) +
          geom_boxplot(width=0.2,aes(fill=Condition)) +scale_y_log10())
