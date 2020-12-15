require(tidyverse)
process_APA <- function(file){
    seqpos <- seq(-10,10,by=1) %>% as.character()
    my.dat <- file %>% read_csv(col_names = F) %>%
        mutate_if(is.character,str_remove,"\\[|\\]") %>%
        mutate_if(is.character,as.numeric) %>%
        as.matrix()
    
    colnames(my.dat) <- seqpos
    rownames(my.dat) <- rev( seqpos )
    my.dat <- my.dat %>% reshape2::melt() %>%
        mutate(file = file) %>%
        mutate(DSB = str_extract(file,"174clived|80best")) %>%
        mutate(Type = str_extract(file,"anchorLeftAsiSI|damaged|notdamaged")) %>%
        mutate(Condition = str_extract(file,"manipA_OHT|manipA_DIvA|manipB_OHT|manipB_DIvA"))
    return(my.dat)
}

#Deux types de matrices : anchorLeft et Damaged/Undamaged
##ANCHORLEFT
files <- c(
    "/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipA/174clived/loops_DIvA_OHT_anchorLeftAsiSI_manipA_OHT/10000/gw/APA.txt",
    "/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipB/174clived/loops_DIvA_OHT_anchorLeftAsiSI_manipB_OHT/10000/gw/APA.txt",
    "/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipA/174clived/loops_DIvA_OHT_anchorLeftAsiSI_manipA_DIvA/10000/gw/APA.txt",
    "/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipB/174clived/loops_DIvA_OHT_anchorLeftAsiSI_manipB_DIvA/10000/gw/APA.txt"
) %>% map(process_APA) %>% bind_rows()

p1 <- files %>% mutate(value = ifelse(value > 1800,1800,value))%>% ggplot(aes(Var2,Var1,fill=value)) + geom_tile() +
    scale_fill_gradient2(low = "white",high = "#EA2027")  +
    theme_classic(base_size = 24) +
    facet_wrap(~Condition,scales="free") +
    ylab("") + xlab("")

pdf("loops_DIvA_OHT_anchorLeftAsiSI_manipA_manipB.pdf",height=18,width=18)
print(p1)
dev.off()
#Merged manip mean
files <- files %>% separate("Condition",into = c("Manip","Condition"))
files <- files %>% group_by(Var1,Var2,Condition) %>% summarise(meanManip = mean(value))
p1 <- files %>% mutate(meanManip = ifelse(meanManip > 2200,2200,meanManip)) %>% ggplot(aes(Var2,Var1,fill=meanManip)) + geom_tile() +
    scale_fill_gradient2(low = "white",high = "#EA2027")  +
    theme_classic(base_size = 24) +
    facet_wrap(~Condition,scales="free") +
    ylab("") + xlab("")

pdf("loops_DIvA_OHT_anchorLeftAsiSI_mean_manipA.pdf",height=9,width=18)
print(p1)
dev.off()

enrichment <- files %>% mutate(pixelLocation = case_when(
    Var1 %in% 6:10 & Var2 %in% -(10:6) ~ "Top-Left",
    # Var1 %in% 9:13 & Var2 %in% 9:13 ~ "Central-Pixel",
    Var1 == 0 & Var2 ==0  ~ "Central-Pixel",
    Var1 %in% -(6:10) & Var2 %in% -(10:6) ~ "Lower-Left",
    TRUE ~ "Others"
)) %>% filter(pixelLocation != "Others") %>% group_by(pixelLocation,Condition) %>% summarise(meanValPixel = mean(meanManip)) %>% 
    spread(key = pixelLocation,value = meanValPixel) %>% 
    mutate(P2ULenrichment = `Central-Pixel`/`Top-Left`) %>%
    mutate(P2LLenrichment = `Central-Pixel`/`Lower-Left`) 

##TADS DAMAGED/UNDAMAGED
my_combi <- expand.grid(c("OHT","DIvA"),c("damaged","notdamaged"),c("174clived","80best"),c("manipA","manipB"))
files2 <- apply(my_combi,1,function(i){
    str_c(
        "/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/",i[4],"/",i[3],
        "/loops_DIvA_OHT_",i[2],"_",i[4],"_",i[1],
        "/10000/gw/APA.txt"
    )
}) %>% map(process_APA) %>% bind_rows()




my_combi <- expand.grid(c("damaged","notdamaged"),c("174clived","80best"),c("manipA","manipB"))

for(j in 1:nrow(my_combi)){
    i <- my_combi[j,]
    mysub <- files2 %>% filter(DSB == i[[2]],Type == i[[1]],str_detect(Condition,as.vector(i[[3]])))
    myquant <- quantile(mysub$value,0.995)
    mysub <- mysub %>% mutate(value = ifelse(value > myquant,myquant,value))
    myp <- mysub %>%
        ggplot(aes(Var2,Var1,fill=value)) + geom_tile() +
        scale_fill_gradient2(low = "white",high = "#EA2027")  +
        theme_classic(base_size = 24) +
        facet_wrap(~Condition,scales="free") +
        ylab("") + xlab("")
    pdf(str_c("loops_DIvA_OHT_",i[[1]],"_",i[[2]],"_",i[[3]],".pdf",sep=""),height=8,width=18)
    print(myp)
    dev.off()
}
#With mean between A and B
my_combi <- expand.grid(c("damaged","notdamaged"),c("174clived","80best"))

for(j in 1:nrow(my_combi)){
    i <- my_combi[j,]
    mysub <- files2 %>% filter(DSB == i[[2]],Type == i[[1]])
    myquant <- quantile(mysub$value,0.995)
    mysub <- mysub %>% mutate(value = ifelse(value > myquant,myquant,value))
    mysub <- mysub %>% dplyr::select(Var1,Var2,Condition,value) %>% separate("Condition",into = c("Manip","Condition"))
    mysub <- mysub %>% group_by(Var1,Var2,Condition) %>% summarise(meanManip = mean(value))
    myp <- mysub %>%
        ggplot(aes(Var2,Var1,fill=meanManip)) + geom_tile() +
        scale_fill_gradient2(low = "white",high = "#EA2027")  +
        theme_classic(base_size = 24) +
        facet_wrap(~Condition,scales="free") +
        ylab("") + xlab("")
    pdf(str_c("loops_DIvA_OHT_",i[[1]],"_",i[[2]],"_meanManip.pdf",sep=""),height=8,width=18)
    print(myp)
    dev.off()
}

#APA CTCF
files <- list.files("APA_CTCF",full.names = T) %>% map(process_APA) %>% bind_rows()

p1 <- files %>% mutate(value = ifelse(value > 150000,150000,value)) %>% ggplot(aes(Var2,Var1,fill=value)) + geom_tile() +
    scale_fill_gradient2(low = "white",high = "#EA2027")  +
    theme_classic(base_size = 24) +
    facet_wrap(~Condition,scales="free") +
    ylab("") + xlab("")

pdf("loops_DIvA_OHT_ctcf_manipA_manipB.pdf",height=18,width=18)
print(p1)
dev.off()

#APA siRNA
#500kb
files <- c(
    "/home/rochevin/Documents/PROJET_INGE/APA/APA_loop_leftsideAsiSI/loops_DIvA_OHT_500kbfromAsiSI5_manipA_DIvA/10000/gw/APA.txt",
    "/home/rochevin/Documents/PROJET_INGE/APA/APA_loop_leftsideAsiSI/loops_DIvA_OHT_500kbfromAsiSI5_manipA_OHT/10000/gw/APA.txt",
    "/home/rochevin/Documents/PROJET_INGE/APA/APA_loop_leftsideAsiSI/loops_DIvA_OHT_500kbfromAsiSI5_manipB_DIvA/10000/gw/APA.txt",
    "/home/rochevin/Documents/PROJET_INGE/APA/APA_loop_leftsideAsiSI/loops_DIvA_OHT_500kbfromAsiSI5_manipB_OHT/10000/gw/APA.txt",
    "/home/rochevin/Documents/PROJET_INGE/APA/APA_loop_leftsideAsiSI/loops_DIvA_OHT_anchorLeftAsiSI_manipA_DIvA/10000/gw/APA.txt",
    "/home/rochevin/Documents/PROJET_INGE/APA/APA_loop_leftsideAsiSI/loops_DIvA_OHT_anchorLeftAsiSI_manipA_OHT/10000/gw/APA.txt",
    "/home/rochevin/Documents/PROJET_INGE/APA/APA_loop_leftsideAsiSI/loops_DIvA_OHT_anchorLeftAsiSI_manipB_DIvA/10000/gw/APA.txt",
    "/home/rochevin/Documents/PROJET_INGE/APA/APA_loop_leftsideAsiSI/loops_DIvA_OHT_anchorLeftAsiSI_manipB_OHT/10000/gw/APA.txt"
) %>% map(process_APA) %>% bind_rows() %>% 
    mutate(Type = str_extract(file,"500kbfromAsiSI|anchorLeftAsiSI"))

p1 <- files %>% filter(Type == "500kbfromAsiSI")%>% mutate(value = ifelse(value > 10000,10000,value))%>% ggplot(aes(Var2,Var1,fill=value)) + geom_tile() +
    scale_fill_gradient2(low = "white",high = "#EA2027")  +
    theme_classic(base_size = 24) +
    facet_wrap(Type~Condition,scales="free") +
    ylab("") + xlab("")

p2 <- files %>% filter(Type == "anchorLeftAsiSI")%>% mutate(value = ifelse(value > 2000,2000,value))%>% ggplot(aes(Var2,Var1,fill=value)) + geom_tile() +
    scale_fill_gradient2(low = "white",high = "#EA2027")  +
    theme_classic(base_size = 24) +
    facet_wrap(Type~Condition,scales="free") +
    ylab("") + xlab("")



pdf("loops_DIvA_OHT_APA_loop_leftsideAsiSI_500kbfromAsiSI_manipA_manipB.pdf",height=36,width=18)
print(cowplot::plot_grid(p1,p2,ncol=1))
dev.off()
