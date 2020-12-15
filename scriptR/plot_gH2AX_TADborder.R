# Raphael Mourad
# 07/02/2019


# Script to make profiles of gH2AX data across TAD borders.



# WORKING DIRECTORY ----------------------------------------------------

setwd("/media/mourad/diskSave/MCF_Toulouse/recherche/LegubeTeam/")

# LIBRARIES ------------------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(HiTC)
library(seqplots)

# Genome 
seqInfohg19=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
seqlen=seqlengths(seqInfohg19)

# Parameters
binSize=1000
win=1000000L # +/- win 


# View points. CHECKED
VP.GR=import.bed("data/4Cseq/4Cviewpoints.bed")
VP.GR=VP.GR[-4]

# DSB. CHECKED
AsiSI="174clived" # "174clived" "HR" "NHEJ"
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
Rad21.GR=import.bw("data/ChIPseq/RAD21/DiVA/HK27MBGX7_RAD21_18s001897-1-1_Clouaire_lane118s001897_sequence_normalized_fixed_step.bw")
Rad21.cov=coverage(Rad21.GR,weight="score")


# Load TAD borders. CHECKED
thres=-0.05
TADborder.GR=import.bed(paste0("data/HiC_Legube/domain/TopDom/border/TADborder_topdom_DIvA_manipA_hg19_50kb_ws10_strong_t",thres,".bed"))
start(TADborder.GR)=end(TADborder.GR)

# TAD border close to VP. CHECKED
TADborderVP.GR=TADborder.GR[nearest(VP.GR,TADborder.GR)]

# TAD border close to DSB. CHECKED
nearestTADborderDSB=nearest(AsiSI.GR,TADborder.GR)
TADborderDSB.GR=TADborder.GR[nearestTADborderDSB[!is.na(nearestTADborderDSB)]]
AsiSI.GR=AsiSI.GR[!is.na(nearestTADborderDSB)] # one AsiSI sites does not have a TAD close


# Average profile on DSB. CHECKED
gH2AX_averageProf=NULL
for(j in 1:length(AsiSI.GR)){
 AsiSIj.GR=AsiSI.GR[j]
 chr=as.character(seqnames(AsiSIj.GR))
 binj.GR=GRanges(chr,IRanges(breakInChunks(totalsize=seqlen[names(seqlen)==chr],chunksize=binSize)))

 TADborderDSBwinj.GR=resize(TADborderDSB.GR[j],width=win*2,fix="center")
 binj_TADborderDSB.GR=binj.GR[subjectHits(findOverlaps(TADborderDSBwinj.GR,binj.GR))]

 distDSBTAD=start(AsiSIj.GR)-start(TADborderDSB.GR[j])

 gH2AXj=viewMeans(Views(gH2AX.cov[[chr]],ranges(binj_TADborderDSB.GR)))
 gH2AXj=lowess(gH2AXj,f=0.02)$y
 if(distDSBTAD<0){
  gH2AXj=gH2AXj[length(gH2AXj):1]
 }

 Rad21j=viewMeans(Views(Rad21.cov[[chr]],ranges(binj_TADborderDSB.GR)))
 Rad21j=lowess(Rad21j,f=0.02)$y
 if(distDSBTAD<0){
  Rad21j=Rad21j[length(Rad21j):1]
 }

 if(j==1){
  gH2AX_averageProf=gH2AXj
  Rad21_averageProf=Rad21j
 }
 if(length(gH2AXj)==length(gH2AX_averageProf)){
  gH2AX_averageProf=gH2AX_averageProf+gH2AXj
  Rad21_averageProf=Rad21_averageProf+Rad21j  
 }
 #plot(gH2AXj)
 print(j)
}


ax=(1:length(binj_TADborderDSB.GR)-(length(binj_TADborderDSB.GR)/2)+0.5)*binSize/1e3
filej=paste0("results/HiC/plotgH2AX/averageProf_gH2AX_TADborders_",AsiSI,".pdf")
pdf(filej,8,8)
par(mfrow = c(2,1))
plot(ax,gH2AX_averageProf,type="l",ylab="gH2AX",xlab="Distance from TAD border (kb)",main="DSB on the right side")
abline(v=0,lty=2,col="red",lwd=2.5)
plot(ax,Rad21_averageProf,type="l",ylab="Rad21",xlab="Distance from TAD border (kb)",main="DSB on the right side")
dev.off()



