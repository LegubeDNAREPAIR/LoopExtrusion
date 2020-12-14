require(tidyverse)
require(plyranges)
bless80 = read_bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")

#Refaire plot 4D avec manip A et B (et plus tard si CTRL)
files <- c(
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipA/80best/loops_verysmall_DIvA_OHT_notdamaged_manipA_DIvA/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipA/notdamaged/","DIvA",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipA/80best/loops_verysmall_DIvA_OHT_notdamaged_manipA_OHT/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipA/notdamaged/","OHT",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipA/80best/loops_verysmall_DIvA_OHT_damaged_manipA_DIvA/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipA/damaged/","DIvA",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipA/80best/loops_verysmall_DIvA_OHT_damaged_manipA_OHT/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipA/damaged/","OHT",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipB/80best/loops_verysmall_DIvA_OHT_notdamaged_manipB_DIvA/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipB/notdamaged/","DIvA",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipB/80best/loops_verysmall_DIvA_OHT_notdamaged_manipB_OHT/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipB/notdamaged/","OHT",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipB/80best/loops_verysmall_DIvA_OHT_damaged_manipB_DIvA/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipB/damaged/","DIvA",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipB/80best/loops_verysmall_DIvA_OHT_damaged_manipB_OHT/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipB/damaged/","OHT",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+")))
) %>%
    map(read_delim," ",col_names = F) %>%
    map(gather) %>%
    bind_rows(.id = "Type") %>%
    separate(Type,into=c("Manip","Type","Condition","seqnames")) %>% 
    mutate(seqnames = str_c("chr",str_extract(seqnames,"[A-Z0-9]+"))) 
files_ratio <- files %>% spread(key= Condition,value=value) %>%
    mutate(ratio = log2(OHT/DIvA)) 
p1 <- files_ratio %>% filter(Manip =="manipA")  %>% mutate(Type = fct_relevel(Type,"notdamaged","damaged")) %>% ggplot(aes(x=Type,y=ratio)) +geom_boxplot() + geom_hline(yintercept = 0,col="red",linetype="dashed")+ facet_wrap(~Manip,ncol=1)
p2 <- files_ratio %>% filter(Manip =="manipB")  %>% mutate(Type = fct_relevel(Type,"notdamaged","damaged"))%>% ggplot(aes(x=Type,y=ratio)) +geom_boxplot() + geom_hline(yintercept = 0,col="red",linetype="dashed")+ facet_wrap(~Manip,ncol=1)

