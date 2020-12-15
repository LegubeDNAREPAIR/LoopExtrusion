# Raphaël MOURAD
# 19/02/2019



# 



# WORKING DIRECTORY ----------------------------------------------------

setwd("/media/mourad/diskSave/MCF_Toulouse/recherche/LegubeTeam")

# LIBRARIES ------------------------------------------------------------



library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(GenomicRanges)
library(Matrix)
library(glmnet)
library(MASS)
library(pryr)
library(HiTC)
library(data.table)
library(ggplot2)
library(biomaRt)
library(circlize)
library(gplots)
library(doMC)
registerDoMC(cores=4)


# PARAMETERS ------------------------------------------------------------

# Genomes
Chr.V=paste0("chr",c(1:22,"X"))
seqInfohg19=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
seqInfohg19=seqInfohg19[Chr.V]
seqlen=seqlengths(seqInfohg19)

# Parameters
binSize=25e3
binSizeTxt=paste0(binSize/1e3,"kb")
distMax=500e3
type.measure="mse"
lambdaVec=exp(-seq(-5,4,0.1))
repet="manipB"


# COMPUTE DIFFERENTIAL PROFILES  ------------------------------------------------------------

# AsiSI
AsiSI="ASIsites_hg19_174_clived_IQ1.5" # "ASIsites_hg19_174_clived_IQ1.5" "BLESS_80best_JunFragPE_Rmdups_pm500bp"
AsiSI.GR=import.bed(paste0("data/AsiSI/",AsiSI,".bed"))
AsiSI.GR=resize(AsiSI.GR,fix="center",width=2e6)


# Load insulation scores
insu_DIvA.GR=import.bedGraph(paste0("data/HiC_Legube/matrix/dump_observed_KR_DIvA_",repet,"_allchr_25kb.is400001.ids50001.insulation.bedGraph"))
insu_OHT.GR=import.bedGraph(paste0("data/HiC_Legube/matrix/dump_observed_KR_OHT_",repet,"_allchr_25kb.is400001.ids50001.insulation.bedGraph"))

# Load TADs
TAD_DIvA.GR=import.bed(paste0("data/HiC_Legube/domain/TopDom/TAD_topdom_DIvA_",repet,"_hg19_50kb_ws5.bed"))
TAD_OHT.GR=import.bed(paste0("data/HiC_Legube/domain/TopDom/TAD_topdom_OHT_",repet,"_hg19_50kb_ws5.bed"))
TAD_DIvA.GRb=unique(GRanges(c(seqnames(TAD_DIvA.GR),seqnames(TAD_DIvA.GR)),
	IRanges(c(start(TAD_DIvA.GR),end(TAD_DIvA.GR)+1),c(start(TAD_DIvA.GR),end(TAD_DIvA.GR)+1))))
TAD_OHT.GRb=unique(GRanges(c(seqnames(TAD_OHT.GR),seqnames(TAD_OHT.GR)),
	IRanges(c(start(TAD_OHT.GR),end(TAD_OHT.GR)+1),c(start(TAD_OHT.GR),end(TAD_OHT.GR)+1))))

TAD_DIvA.GRb_d=TAD_DIvA.GRb[unique(queryHits(findOverlaps(TAD_DIvA.GRb,AsiSI.GR)))]
TAD_OHT.GRb_d=TAD_OHT.GRb[unique(queryHits(findOverlaps(TAD_OHT.GRb,AsiSI.GR)))]
TAD_DIvA.GRb_nd=TAD_DIvA.GRb[-unique(queryHits(findOverlaps(TAD_DIvA.GRb,AsiSI.GR)))]
TAD_OHT.GRb_nd=TAD_OHT.GRb[-unique(queryHits(findOverlaps(TAD_OHT.GRb,AsiSI.GR)))]

# Profile
win=10
#olid=findOverlaps(TAD_DIvA.GRb_d,insu_DIvA.GR)
olid=findOverlaps(unique(c(TAD_DIvA.GRb_d,TAD_OHT.GRb_d)),insu_DIvA.GR)
matid=sapply(subjectHits(olid),function(x){(x-win):(x+win)})
profDIvA_d=t(apply(matid,2,function(x){insu_DIvA.GR$score[x]}))
profOHT_d=t(apply(matid,2,function(x){insu_OHT.GR$score[x]}))
profRatio_d=profOHT_d-profDIvA_d

#olind=findOverlaps(TAD_DIvA.GRb_nd,insu_DIvA.GR)
olind=findOverlaps(unique(c(TAD_DIvA.GRb_nd,TAD_OHT.GRb_nd)),insu_DIvA.GR)
matind=sapply(subjectHits(olind),function(x){(x-win):(x+win)})
profDIvA_nd=t(apply(matind,2,function(x){insu_DIvA.GR$score[x]}))
profOHT_nd=t(apply(matind,2,function(x){insu_OHT.GR$score[x]}))
profRatio_nd=profOHT_nd-profDIvA_nd

col=colorRampPalette(c("#2980b9","black","#f1c40f"))(255)
breaks=seq(-0.1, 0.1, length.out=256)
distPlot=paste0(seq(-(binSize/1e3*10),(binSize/1e3*10),binSize/1e3),"kb")
dist=seq(-(binSize/1e3*10),(binSize/1e3*10),binSize/1e3)

# Average prof
wt_d=wilcox.test(profDIvA_d[,win+1],profOHT_d[,win+1])
#FC_d=median(profOHT_d[,win+1])/median(profDIvA_d[,win+1])
FC_d=median(profOHT_d[,win+1])-median(profDIvA_d[,win+1])
boxplot(profDIvA_d[,win+1],profOHT_d[,win+1],ylim=c(-1,1))
file_aveProf_insu_d="results/HiC/insuProfile/averageProf_insu_OHT_VS_DIvA_damaged.pdf"
pdf(file_aveProf_insu_d,5,5)
plot(dist,apply(profDIvA_d,2,median),xlab="Distance to TAD border (kb)",ylab="Insulation score (OHT-DIvA)",col="blue",type="l",ylim=c(-0.4,.2))
lines(dist,apply(profOHT_d,2,median),col="red",type="l")
legend("topright",legend=c("DIvA","OHT"),col=c("blue","red"),lty=1)
dev.off()

wt_nd=wilcox.test(profDIvA_nd[,win+1],profOHT_nd[,win+1])
FC_nd=median(profOHT_nd[,win+1])-median(profDIvA_nd[,win+1])
boxplot(profDIvA_nd[,win+1],profOHT_nd[,win+1],ylim=c(-1,1))
file_aveProf_insu_nd="results/HiC/insuProfile/averageProf_insu_OHT_VS_DIvA_notdamaged.pdf"
pdf(file_aveProf_insu_nd,5,5)
plot(dist,apply(profDIvA_nd,2,median),xlab="Distance to TAD border (kb)",ylab="Insulation score (OHT-DIvA)",col="blue",type="l",ylim=c(-0.4,.2))
lines(dist,apply(profOHT_nd,2,median),col="red",type="l")
legend("topright",legend=c("DIvA","OHT"),col=c("blue","red"),lty=1)
dev.off()



##################
# Insulation score 10kb at exact DSB site
# Load insulation scores
insu_DIvA.GR=import.bedGraph(paste0("data/HiC_Legube/matrix/insu10kb/dump_observed_KR_DIvA_",repet,"_allchr_10kb.is160001.ids20001.insulation.bedGraph"))
insu_OHT.GR=import.bedGraph(paste0("data/HiC_Legube/matrix/insu10kb/dump_observed_KR_OHT_",repet,"_allchr_10kb.is160001.ids20001.insulation.bedGraph"))

# CTCF
dataPeak=read.table("data/ChIPseq/CTCF/CTCFalfDIvA_vs_INPUT_peaks.narrowPeak",sep="\t",header=F)
CTCFpeak.GR=GRanges(dataPeak[,1],IRanges(dataPeak[,2],dataPeak[,3]))
CTCFpeak.GR=CTCFpeak.GR[start(CTCFpeak.GR)>1e6]
CTCFpeak_d.GR=subsetByOverlaps(CTCFpeak.GR,resize(AsiSI.GR,width=1e6, fix="center"))

center.GR=anchorDIvA.GR # AsiSI.GR or CTCFpeak.GR or CTCFpeak_d.GR anchorDIvA1.GR

# Profile
binSize=10e3
win=100
olid=findOverlaps(center.GR,insu_DIvA.GR)
matid=sapply(unique(subjectHits(olid)),function(x){(x-win):(x+win)})
profDIvA_d=t(apply(matid,2,function(x){insu_DIvA.GR$score[x]}))
profOHT_d=t(apply(matid,2,function(x){insu_OHT.GR$score[x]}))
profDIvA_d=na.omit(profDIvA_d)
profOHT_d=na.omit(profOHT_d)

col=colorRampPalette(c("#2980b9","black","#f1c40f"))(255)
breaks=seq(-0.1, 0.1, length.out=256)
dist=seq(-(binSize/1e3*win),(binSize/1e3*win),binSize/1e3)

# Average prof
wt_d=wilcox.test(profDIvA_d[,win+1],profOHT_d[,win+1])
FC_d=median(profOHT_d[,win+1])-median(profDIvA_d[,win+1])
#boxplot(profDIvA_d[,win+1],profOHT_d[,win+1],ylim=c(-1,1))
file_aveProf_insu_d="results/HiC/insuProfile/averageProf_insu_OHT_VS_DIvA_atAsiSI.pdf"
pdf(file_aveProf_insu_d,5,5)
plot(dist,apply(profDIvA_d,2,median),xlab="Distance to center (kb)",ylab="Insulation score (OHT-DIvA)",col="blue",type="l",ylim=c(0,.2))
lines(dist,apply(profOHT_d,2,median),col="red",type="l")
legend("topright",legend=c("DIvA","OHT"),col=c("blue","red"),lty=1)
dev.off()


##################
# Insulation score 10kb at loop anchors depending on loop size
# loop anchors

dataLoopDIvA=read.table(paste0("data/HiC_Legube/loop/hiccups/DIvA_",repet,"/merged_loops"),sep="\t",header=T)
dataLoopDIvA$size=(dataLoopDIvA[,6]-dataLoopDIvA[,3])
sizeMin=c(seq(100e3,1.9e6,by=100e3),seq(3e6,9e6,by=1e6))
sizeMax=c(seq(200e3,2e6,by=100e3),seq(4e6,10e6,by=1e6))
mean_insu_DIvA=NULL
mean_insu_OHT=NULL
for(i in 1:length(sizeMin)){
 dataLoopDIvAi=dataLoopDIvA[dataLoopDIvA$size>=sizeMin[i] & dataLoopDIvA$size<sizeMax[i],]
 anchorDIvA.GRi=GRanges(c(paste0("chr",dataLoopDIvAi[,1]),paste0("chr",dataLoopDIvAi[,1])),IRanges(c((dataLoopDIvAi[,2]+dataLoopDIvAi[,3])/2,(dataLoopDIvAi[,5]+dataLoopDIvAi[,6])/2),c((dataLoopDIvAi[,2]+dataLoopDIvAi[,3])/2,(dataLoopDIvAi[,5]+dataLoopDIvAi[,6])/2)))

 olid=findOverlaps(anchorDIvA.GRi,insu_DIvA.GR)
 mean_insu_DIvA=c(mean_insu_DIvA,mean(insu_DIvA.GR$score[subjectHits(olid)]))
 mean_insu_OHT=c(mean_insu_OHT,mean(insu_OHT.GR$score[subjectHits(olid)]))
}

file_insu_loopsize="results/HiC/insuProfile/insu_OHT_VS_DIvA_loopsize.pdf"
pdf(file_insu_loopsize)
plot((sizeMin+sizeMax)/2e3,mean_insu_DIvA,type="l",xlab="Loop size (kb)",ylab="Insulation score",col="blue",
		ylim=c(-0.6,0.1), main="Insulation score at loop anchors depending on loop size")
lines((sizeMin+sizeMax)/2e3,mean_insu_OHT,type="l",xlab="Loop size (kb)",ylab="Insulation score",col="red")
legend("topright",legend=c("DIvA","OHT"),col=c("blue","red"),lty=1)
dev.off()







