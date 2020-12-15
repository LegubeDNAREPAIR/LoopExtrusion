# Raphael Mourad
# 07/02/2019


# Script to plot 4C data with other data.



# WORKING DIRECTORY ----------------------------------------------------

setwd("/media/mourad/diskSave/MCF_Toulouse/recherche/LegubeTeam/")

# LIBRARIES ------------------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(HiTC)
library(scales)

# Genome 
seqInfohg19=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
seqlen=seqlengths(seqInfohg19)

# Parameters
win=2e6
reskb="50kb"
res=50000

# View points. CHECKED
VP.GR=import.bed("data/4Cseq/4Cviewpoints.bed")
VP.GR=VP.GR[-4]

# Load ChIP-seq. CHECKED
gH2AX.GR=import.wig("../../data/BIGWIG/GAM.clean_hg19.wig")
gH2AX.cov=coverage(gH2AX.GR,weight="score")
Rad21.GR=import.bw("../../data/BIGWIG/HLNKYBGXC_Pool_ChIP-seq_legube_19s004478-1-1_Clouaire_lane1Rad21DIvA_sequence.exp_spikeinfactor.bw")
Rad21.cov=coverage(Rad21.GR,weight="score")
CTCF.GR=import.bw("../../data/BIGWIG/H5HJKBGXC_ChIP-seq_thomas_19s003259-1-1_Clouaire_lane1CTCFalfmin_sequence_normalized.bw")
CTCF.cov=coverage(CTCF.GR,weight="score")
dataPeak=read.table("../../data/peak_calling/CTCF/MACS/CTCFalfDIvA_vs_INPUT_peaks.narrowPeak",sep="\t",header=F)
dataPeak=dataPeak[dataPeak[,7]>quantile(dataPeak[,7],0),]
CTCFpeak.GR=GRanges(dataPeak[,1],IRanges(dataPeak[,2],dataPeak[,3]))

# Load CTCF motifs hg19
dataCTCFmotif=read.table("../../data/ctcf_motifs.gff",sep="\t")
CTCFmotif.GR=GRanges(dataCTCFmotif[,1],IRanges(dataCTCFmotif[,4],dataCTCFmotif[,5]),score=dataCTCFmotif[,6],strand=dataCTCFmotif[,7])
CTCFmotif.GR=CTCFmotif.GR[CTCFmotif.GR$score>quantile(CTCFmotif.GR$score,0.4)]
CTCFmotifp.cov=coverage(CTCFmotif.GR[strand(CTCFmotif.GR)=='+'])
CTCFmotifm.cov=coverage(CTCFmotif.GR[strand(CTCFmotif.GR)=='-'])

# Peak orientation
CTCFmotifp.GR=CTCFmotif.GR[strand(CTCFmotif.GR)=='+']
CTCFmotifm.GR=CTCFmotif.GR[strand(CTCFmotif.GR)=='-']
CTCFpeakp.GR=subsetByOverlaps(CTCFpeak.GR,CTCFmotifp.GR)
CTCFpeakm.GR=subsetByOverlaps(CTCFpeak.GR,CTCFmotifm.GR)

# Load 4Cseq. Same order as view points. CHECKED
FourCfiles=c("segToFrag_minusOHT_ctrl_rich_smoothed_11FragsPerWin.bedGraph.gz",
	"segToFrag_minusOHT_chr17_AsiSI_smoothed_11FragsPerWin.bedGraph.gz",
	"segToFrag_minusOHT_chr1_AsiSI_smoothed_11FragsPerWin.bedGraph.gz",
	"segToFrag_minusOHT_chr20_AsiSI_smoothed_11FragsPerWin.bedGraph.gz",
	"segToFrag_minusOHT_chr21_AsiSI_smoothed_11FragsPerWin.bedGraph.gz")

# Load TAD borders. CHECKED
thres=-0.05
TADborder.GR=import.bed(paste0("data/HiC_Legube/domain/TopDom/border/TADborder_topdom_DIvA_manipA_hg19_50kb_ws10_strong_t",thres,".bed"))
start(TADborder.GR)=end(TADborder.GR)

# Process Hi-C. CHECKED
load(paste0("data/HiC_Legube/HTC/HTC_observed_KR_DIvA_manipA_intrachr_",reskb,".RData"))
HTCjsub_list=list()
for(j in 1:length(VP.GR)){
 VPj.GR=VP.GR[j]
 VPwinj.GR=resize(VP.GR[j],width=win*2,fix="center")
 chr=as.character(seqnames(VPj.GR))
 HTCjsub=extractRegion(HTCL[[paste0(chr,chr)]], c(1,2), chr=chr, from=start(VPwinj.GR), to=end(VPwinj.GR))
 HTCjsub_list[[j]]=HTCjsub
 print(j)
}


# Load insulation score
insu.GR=import(paste0("data/HiC_Legube/matrix/insu",reskb,"/dump_observed_KR_DIvA_manipA_allchr_",reskb,".is800001.ids100001.insulation.bedGraph"))

# Plot 
for(j in 1:length(VP.GR)){

 VPj.GR=VP.GR[j]
 VPwinj.GR=resize(VP.GR[j],width=win*2,fix="center")
 chr=as.character(seqnames(VPj.GR))
 binj.GR=GRanges(chr,IRanges(breakInChunks(totalsize=seqlen[names(seqlen)==chr],chunksize=res)))
 binj_VP.GR=binj.GR[subjectHits(findOverlaps(VPwinj.GR,binj.GR))]
 bin1kj.GR=GRanges(chr,IRanges(breakInChunks(totalsize=seqlen[names(seqlen)==chr],chunksize=50000)))
 bin1kj_VP.GR=bin1kj.GR[subjectHits(findOverlaps(VPwinj.GR,bin1kj.GR))]
 bin10j.GR=GRanges(chr,IRanges(breakInChunks(totalsize=seqlen[names(seqlen)==chr],chunksize=10)))
 bin10j_VP.GR=bin10j.GR[subjectHits(findOverlaps(VPwinj.GR,bin10j.GR))]

 # Load insulation score. GM12878 used
 insu.GRj=insu.GR[subjectHits(findOverlaps(VPwinj.GR,binj.GR))]

 # Load 4Cseq
 FourC.GR=import.bedGraph(paste0("data/4Cseq/minusOHT/",FourCfiles[j]))
 FourC.cov=coverage(FourC.GR,weight="score")
 FourCj=viewMeans(Views(FourC.cov[[chr]],ranges(binj_VP.GR)))
 FourCj[FourCj>quantile(FourCj,0.985)]=quantile(FourCj,0.985) # thresholding of 4C data to compare with gH2AX

 # ChIP-seq
 gH2AXj=viewMeans(Views(gH2AX.cov[[chr]],ranges(binj_VP.GR)))
 Rad21j=viewMeans(Views(Rad21.cov[[chr]],ranges(bin1kj_VP.GR)))
 CTCFj=viewMeans(Views(CTCF.cov[[chr]],ranges(bin1kj_VP.GR)))
 CTCFpj1k=countOverlaps(bin1kj_VP.GR,CTCFpeakp.GR)
 CTCFmj1k=countOverlaps(bin1kj_VP.GR,CTCFpeakm.GR)
 CTCFpj10=countOverlaps(bin10j_VP.GR,CTCFpeakp.GR)
 CTCFmj10=countOverlaps(bin10j_VP.GR,CTCFpeakm.GR)

 # CTCF motifs
 CTCFmotifpj=Views(CTCFmotifp.cov[[chr]],ranges(VPwinj.GR)) 
 CTCFmotifmj=Views(CTCFmotifm.cov[[chr]],ranges(VPwinj.GR)) 

 # TAD borders
 #TADborderj=viewMeans(Views(coverage(TADborder.GR)[[chr]],ranges(binj_VP.GR)))
 #TADborderj[TADborderj>0]=1
 TADborderj=viewMeans(Views(coverage(TADborder.GR)[[chr]],ranges(bin1kj_VP.GR)))
 TADborderj[TADborderj>0]=1

 # Plots
 ax=(1:length(binj_VP.GR)-(length(binj_VP.GR)/2)+0.5)*res/1e3
 ax1k=(1:length(bin1kj_VP.GR)-(length(bin1kj_VP.GR)/2)+0.5)
 ax10=(1:length(bin10j_VP.GR)-(length(bin10j_VP.GR)/2)+0.5)
 filej=paste0("results/4Cseq/plot4Cseq/plot4Cseq_othertracks_",seqnames(VPj.GR),"_pos",start(VPj.GR),".pdf")
 pdf(filej,4,6)
 par(mfrow = c(8,1))
 par(mar=c(1,5,1,1))
 plot(ax,gH2AXj,type='l',ylab="gH2AX",xlab="")
 plot(ax,FourCj,type='l',ylab="4C",xlab="")
 plot(ax,insu.GRj$score,type='l',ylab="Insulation score",xlab="")
 #plot(ax,TADborderj,type='l',ylab="TAD border",xlab="")
 plot(ax1k,rep(0,length(ax1k)),ylim=c(0,1),type='l',ylab="TAD border",xlab="",col="white",yaxt="n")
 abline(v=(which(TADborderj==1)-length(ax1k)/2),col="red")
 plot(ax1k,Rad21j,type='l',ylab="Rad21",xlab="")
 plot(ax1k,CTCFj,type='l',ylab="CTCF",xlab="")
 plot(ax1k,CTCFpj1k,type='l',ylab="CTCF peak +\n density",ylim=c(0,3),xlab="")
 plot(ax1k,CTCFmj1k,type='l',ylab="CTCF peak -\n density",ylim=c(0,3),xlab="")
 plot(ax1k,CTCFpj1k,type='l',ylab="CTCF peak\n density (10kb)",ylim=c(0,3),xlab="",col=alpha("blue",0.5))
 lines(ax1k,CTCFmj1k,type='l',col=alpha("red",0.5))
 plot(ax10/1e3,CTCFpj10,type='l',ylab="CTCF peak\n (no binning)",ylim=c(0,3),xlab="",col=alpha("blue",0.5))
 lines(ax10/1e3,CTCFmj10,type='l',col=alpha("red",0.5))

 #plot((-(length(as.vector(CTCFmotifmj)[[1]])/2):(length(as.vector(CTCFmotifmj)[[1]])/2))[-1]/1000,as.vector(CTCFmotifmj)[[1]],type='l',ylab="CTCF motif +",xlab="")
 #plot((-(length(as.vector(CTCFmotifpj)[[1]])/2):(length(as.vector(CTCFmotifpj)[[1]])/2))[-1]/1000,as.vector(CTCFmotifpj)[[1]],type='l',ylab="CTCF motif -",xlab="")
 dev.off()

 filej=paste0("results/4Cseq/plot4Cseq/plotHiC_othertracks_",seqnames(VPj.GR),"_pos",start(VPj.GR),".pdf")
 pdf(filej,8,4)
 mapC(HTCjsub_list[[j]],log.data=T) # use GM12878
 dev.off()

}





