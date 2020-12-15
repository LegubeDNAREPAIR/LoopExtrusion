require(tidyverse)
time_gamma <- "../data/speed_gamma.csv" %>% read_csv2() %>% mutate(DSB = str_c("DSB",1:dplyr::n()))

p<- time_gamma %>%
    gather(key = Time,value = value,-DSB) %>% 
    mutate(Time = recode(Time, `60` = "3600", `300` = "300",`900`="900")) %>% 
    mutate(Time = factor(Time,levels = c("300","900","3600"))) %>% 
    ggplot(aes(x=Time,y=value,fill=Time)) + geom_boxplot() + theme_classic() + theme(legend.position = "none") +
    ylab("Distance (bp)")+xlab("Time(s)")
print(p)

p2 <- time_gamma %>%
    mutate(`0` = 0) %>% 
    gather(key = Time,value = value,-DSB) %>% 
    mutate(Time = recode(Time, `60` = "3600", `300` = "300",`900`="900")) %>% 
    mutate(Time = as.numeric(Time)) %>% 
    ggplot(aes(x=Time,y=value)) +
    stat_summary(fun.y = mean, geom = "point") + 
    stat_summary(fun.y = mean,geom = "line",group=1) +
    stat_summary(fun.data = mean_se, geom = "errorbar",width=5) +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(expand = c(0, 0),breaks = c(300,900,3600),limits = c(0,3800)) +
    scale_y_continuous(expand = c(0, 0)) +
    ylab("Distance (bp)")+xlab("Time(s)") +
    theme_classic()
print(p2)


#Compute pente
Lm_me <- function(x,by_zero=T){
    if(by_zero)
        my_lm <- summary(lm(value~0+Time,data=x))$coefficient[1]
    else
        my_lm <-  summary(lm(value~Time,data=x))$coefficient[2,1]
    return(my_lm)
}
##Let's go:
ccdata <- time_gamma %>%
    gather(key = Time,value = value,-DSB) %>% 
    mutate(Time = recode(Time, `60` = "3600", `300` = "300",`900`="900")) %>% 
    mutate(Time = as.numeric(Time)) %>% group_by(DSB) %>% nest() %>% mutate(lm_by_zero = map_dbl(data,Lm_me)) %>% mutate(lm_no_zero = map_dbl(data,Lm_me,by_zero=F)) %>% 
    unnest(data) %>% 
    spread(key = Time,value=value)

p7 <- ccdata %>% dplyr::select(DSB,lm_by_zero,lm_no_zero) %>% gather(key = Type,value = Vitesse,-DSB) %>% ggplot(aes(x=Type,y=Vitesse)) + geom_boxplot() + 
    theme_classic()
print(p7)
