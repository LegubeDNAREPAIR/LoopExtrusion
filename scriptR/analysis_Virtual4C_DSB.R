# Raphael Mourad
# 28/08/2019


# GOOD RESULTS!


# WORKING DIRECTORY ----------------------------------------------------

setwd("/media/mourad/diskSave/MCF_Toulouse/recherche/LegubeTeam/")
library(Matrix)
library(MASS)
library(HiTC)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(HiCblock)
library(gplots)
library(metap)


# LOAD DATA AND RESULTS ----------------------------------------------------------

# Exp
expe="manipA"
obserOE="observed" # "observed" "OE"
AsiSI="80best_shiftleft200kb" # "174clived" "HR" "NHEJ"

# Resolution
binSize=100e3
binSizeTxt=paste0(binSize/1e3,"kb")
win=1e6
winbin=win/binSize

# Chromosomes
SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
Chr.V=paste0("chr",c(1:22,"X"))
seqlen=seqlengths(SeqInfo[Chr.V])

# AsiSI
if(AsiSI=="174clived"){
 AsiSI.GR=import.bed("data/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed")
}else if(AsiSI=="HR"){
 AsiSI.GR=import.bed("data/AsiSI/BLESS_HR_JunFragPE_Rmdups_pm500bp.bed")
}else if(AsiSI=="NHEJ"){
 AsiSI.GR=import.bed("data/AsiSI/BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed")
}else if(AsiSI=="80random"){
 AsiSI.GR=import.bed("data/AsiSI/80random.bed")
}else if(AsiSI=="80best"){
 AsiSI.GR=import.bed("data/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
}else if(AsiSI=="80best_shiftright200kb"){
 AsiSI.GR=import.bed("data/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp_shiftright200kb.bed")
}else if(AsiSI=="80best_shiftleft200kb"){
 AsiSI.GR=import.bed("data/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp_shiftleft200kb.bed")
}


# Load ChIP-seq. CHECKED
file_gH2AX="data/ChIPseq/gH2AX/OHT4h/GAM.clean_hg19.wig"
gH2AX.GR=import.wig(file_gH2AX)
gH2AX.cov=coverage(gH2AX.GR,weight="score")


# Load HiTC per chr
fileDIvA=paste0("data/HiC_Legube/HTC/HTC_",obserOE,"_KR_DIvA_",expe,"_intrachr_",binSizeTxt,".RData")
load(fileDIvA)
HTCL_DIvA=HTCL
fileOHT=paste0("data/HiC_Legube/HTC/HTC_",obserOE,"_KR_OHT_",expe,"_intrachr_",binSizeTxt,".RData")
load(fileOHT)
HTCL_OHT=HTCL

# Extract Hi-C data at DSBs
V4C_profDIvA=NULL
V4C_profOHT=NULL
gH2AX_vec=rep(NA,length(AsiSI.GR))
for(j in 1:length(AsiSI.GR)){
 AsiSIj.GR=AsiSI.GR[j]
 chr=as.character(seqnames(AsiSIj.GR))
 HTC_DIvAj=HTCL_DIvA[[paste0(chr,chr)]]
 HTC_OHTj=HTCL_OHT[[paste0(chr,chr)]]
 facNormj=sum(intdata(HTC_DIvAj))/sum(intdata(HTC_OHTj))
 bin.GRj=x_intervals(HTC_DIvAj)
 idxbinj=which(countOverlaps(bin.GRj,AsiSIj.GR)>0)

 if(length(idxbinj)>0){
 gH2AX_vec[j]=viewMeans(Views(gH2AX.cov[[chr]],ranges(resize(AsiSIj.GR,win,fix="center"))))

 if((idxbinj+winbin)>ncol(intdata(HTC_DIvAj))){
  V4C_profDIvAj=intdata(HTC_DIvAj)[idxbinj,(idxbinj-winbin):ncol(intdata(HTC_DIvAj))]
  V4C_profDIvAj=c(V4C_profDIvAj,rep(0,winbin*2+1-length(V4C_profDIvAj)))
 }else if((idxbinj-winbin)<1){
  V4C_profDIvAj=intdata(HTC_DIvAj)[idxbinj,1:(idxbinj+winbin)]
  V4C_profDIvAj=c(rep(0,winbin*2+1-length(V4C_profDIvAj)),V4C_profDIvAj)
 }else{
  V4C_profDIvAj=intdata(HTC_DIvAj)[idxbinj,(idxbinj-winbin):(idxbinj+winbin)]
 }

 if((idxbinj+winbin)>ncol(intdata(HTC_OHTj))){
  V4C_profOHTj=intdata(HTC_OHTj)[idxbinj,(idxbinj-winbin):ncol(intdata(HTC_OHTj))]
  V4C_profOHTj=c(V4C_profOHTj,rep(0,winbin*2+1-length(V4C_profOHTj)))
 }else if((idxbinj-winbin)<1){
  V4C_profOHTj=intdata(HTC_OHTj)[idxbinj,1:(idxbinj+winbin)]
  V4C_profOHTj=c(rep(0,winbin*2+1-length(V4C_profOHTj)),V4C_profOHTj)
 }else{
  V4C_profOHTj=intdata(HTC_OHTj)[idxbinj,(idxbinj-winbin):(idxbinj+winbin)]
 }

 V4C_profDIvA=rbind(V4C_profDIvA,V4C_profDIvAj)
 V4C_profOHT=rbind(V4C_profOHT,round(V4C_profOHTj*facNormj))
 }
}

binPos=seq(-win,win,binSize)/1e3

V4C_profDIvA=V4C_profDIvA[order(gH2AX_vec,decreasing=T),]
V4C_profDIvAlog=log(V4C_profDIvA)
V4C_profDIvAlog[V4C_profDIvAlog==-Inf]=0
file_heatmapDIvA=paste0("results/HiC/V4C/heatmapDIvA_",binSizeTxt,"_",AsiSI,".pdf")
pdf(file_heatmapDIvA,8,16)
heatmap.2(V4C_profDIvAlog,Rowv=F,Colv=F,trace="none",labRow=F,labCol=binPos,xlab="Distance to DSB (kb)")
dev.off()

V4C_profOHT=V4C_profOHT[order(gH2AX_vec,decreasing=T),]
V4C_profOHTlog=log(V4C_profOHT)
V4C_profOHTlog[V4C_profOHTlog==-Inf]=0
file_heatmapOHT=paste0("results/HiC/V4C/heatmapOHT_",binSizeTxt,"_",AsiSI,".pdf")
pdf(file_heatmapOHT,8,16)
heatmap.2(V4C_profOHTlog,Rowv=F,Colv=F,trace="none",labRow=F,labCol=binPos,xlab="Distance to DSB (kb)")
dev.off()

colfunc <- colorRampPalette(c("blue", "white", "red"))
V4C_profDiff=log2((V4C_profOHT+1)/(V4C_profDIvA+1))
file_heatmapDiff=paste0("results/HiC/V4C/heatmapDiff_",binSizeTxt,"_",AsiSI,".pdf")
pdf(file_heatmapDiff,8,16)
heatmap.2(V4C_profDiff,Rowv=F,Colv=F,trace="none",labRow=F,labCol=binPos,
	col=colfunc(51),xlab="Distance to DSB (kb)",
	breaks=seq(quantile(V4C_profDiff,.01), quantile(V4C_profDiff,.99), length.out=52))
dev.off()

# Compute p-values
pval=sapply(1:length(binPos),function(x){wilcox.test(V4C_profDIvA[,x],V4C_profOHT[,x],paired=T)$p.value})
V4C_profMeanDiff=colMeans(V4C_profOHT)/colMeans(V4C_profDIvA)
V4C_profMeanDiff[abs(binPos)>100 & abs(binPos)<500]

pvalShortDist=pval[abs(binPos)>100 & abs(binPos)<500]
pShortDist=sumlog(pvalShortDist)$p

pvalLongDist=pval[abs(binPos)>500 & abs(binPos)<1000]
pLongDist=sumlog(pvalLongDist)$p


file_averageProf=paste0("results/HiC/V4C/averageProf_",binSizeTxt,"_",AsiSI,".pdf")
pdf(file_averageProf,5,7)
par(mfrow=c(3,1))
par(mar=c(4,4,1,1))
plot(binPos,colMeans(V4C_profDIvA),type="l",log="y",xlab="Distance to DSB (kb)",ylab="Hi-C count")
#plot(binPos,colMeans(V4C_profDIvA),type="l",log="y",xlab="Distance to DSB (kb)",ylab="Hi-C count",
#	main=paste0("p(100-500kb)=",format(pShortDist)," p(500-1000kb)=",format(pLongDist)))
lines(binPos,colMeans(V4C_profOHT),col="red")
legend("topright",legend=c("DIvA","OHT"),col=c("black","red"),lty=1)
plot(binPos,colMeans(V4C_profOHT)/colMeans(V4C_profDIvA),type="l",log="y",xlab="Distance to DSB (kb)",ylab="Fold-change (OHT/DIvA)")
abline(h=1,col="blue")
plot(binPos,pval,xlab="Distance to DSB (kb)",ylab="p-value",type="l",log="y")
dev.off()




