library(Matrix)
library(MASS)
library(HiTC)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
#bash code
# expe="DIvA_manipA OHT_manipA DIvA_manipB OHT_manipB"
# obserOE="observed"
# res=100000
# reskb="100kb"
# for exp in $expe; do for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do  java -jar programs/juicer/juicer_tools.jar dump $obserOE KR /mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Hi-C/HiC/hicfiles/${exp}/inter.hic $chr $chr BP $res data/KR/dump_${obserOE}_KR_${exp}_chr${chr}_${reskb}.txt; gzip data/KR/dump_${obserOE}_KR_${exp}_chr${chr}_${reskb}.txt; done; done

# LOAD DATA AND RESULTS ----------------------------------------------------------

# Exp
obserOE="observed" # "observed" "OE"
expe=c("DIvA_manipA","OHT_manipA","DIvA_manipB","OHT_manipB")
# expe=c("HiC_Legube_DIvA_rep1_with_1kb", "HiC_Legube_DIvA_rep2_with_1kb", "HiC_Legube_OHT_rep1_with_1kb", "HiC_Legube_OHT_rep2_with_1kb")

# Resolution
binSize=25000
binSizeTxt=paste0(binSize/1e3,"kb")


# Chromosomes
SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
Chr.V=paste0("chr",c(1:22,"X"))
seqlen=seqlengths(SeqInfo[Chr.V])


# Build HiTC per chr
for(k in 1:length(expe)){
    HTCl=list()
    for(i in 1:length(Chr.V)){
        
        # Load data
        chrendi=seqlengths(SeqInfo[Chr.V[i]])
        binendi=ceiling(chrendi/binSize)
        
        fileCounti=paste0("data/KR/dump_",obserOE,"_KR_",expe[k],"_",Chr.V[i],"_",Chr.V[i],"_",binSizeTxt,".txt.gz")
        datai=as.matrix(fread(paste0("zcat ",fileCounti),sep='\t',header=F))
        if(nrow(datai) == 0)
            next
        datai=datai[!is.na(datai[,3]),]
        datai=rbind(datai,c((binendi-1)*binSize,(binendi-1)*binSize,0))
        if(obserOE=="observed"){
            datai[,3]=round(datai[,3])
        }else{
            datai[,3]=round(datai[,3],2)
        }
        
        data.Mati=sparseMatrix(i=(datai[,1]/binSize)+1,j=(datai[,2]/binSize)+1,x=datai[,3])
        data.MatSymi=data.Mati+t(data.Mati)
        diag(data.MatSymi)=diag(data.MatSymi)/2
        rm(data.Mati,datai)
        
        starti=seq(1,chrendi,by=binSize)
        endi=c(seq(binSize,chrendi,by=binSize),chrendi)
        chri.GR=GRanges(Chr.V[i],IRanges(starti,endi))
        names(chri.GR)=paste0("bin",1:length(chri.GR))
        
        HTC=HTCexp(data.MatSymi,chri.GR,chri.GR)
        HTCl[[i]]=HTC
        fileouti=paste0("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/data/HTC/OE/HTC_",obserOE,"_KR_",expe[k],"_",Chr.V[i],"_",binSizeTxt,".RData")
        save(HTC,file=fileouti)
        
        rm(data.MatSymi,chri.GR)
    }
    #HTCL=HTClist(HTCl)
    # fileouti=paste0("/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/HiC/HTC/HTC_",obserOE,"_KR_",expe[k],"_intrachr_",binSizeTxt,".RData")
    # save(HTCL,file=fileouti)
}









