## need baitCoords + regToExclude. Take the mean scores of fragments falling into the baitCoords +/- extSize around the bait, fragments into regToExclude not taken into account.
## return the normalised scores.

options(stringsAsFactors = FALSE)
options(scipen=999) # disable the scientific notation



#smoothData_adapted.R : adapted from smoothData (https://github.com/bbcf/bbcfutils/blob/master/R/smoothData.R) to be used in a snakemake pipeline

.libPaths(c(.libPaths(),"/home/rochevin/R-packages/"))
setwd(snakemake@params[["dir"]]) ;
require(dplyr)
options(stringsAsFactors = FALSE)
options(scipen=999) # disable the scientific notation



fragsFile <- snakemake@input[[1]]
normFile <- snakemake@output[[1]]
baitCoords <- snakemake@params[["viewpoint"]]
extSize <- as.numeric(snakemake@params[["extSize"]])
curName <- snakemake@params[["curName"]]
print(paste("curName=",curName,sep=""))
#regToExclude="chr2:1000-1000"; #previously default is a fake region which does not concern any fragments
regToExclude <- snakemake@params[["regToExclude"]]

# might be passed as parameters
nskip=1
curSep="\t"

chr_regToInclude=unlist(strsplit(baitCoords,split=":"))[1]
coord_baitCoords=as.numeric(unlist(strsplit(unlist(strsplit(baitCoords,split=":"))[2],"-")))
mid_baitCoords=round(sum(coord_baitCoords)/2)
if(mid_baitCoords-extSize>0){regStart=mid_baitCoords-extSize}else{regStart=coord_baitCoords[1]} ## make sure the extended region has a positive start coord
coord_regToInclude=c(regStart,mid_baitCoords+extSize)
names(coord_regToInclude)=c("start","end")

## #default is set to bait +/-5kb
if(regToExclude=="NA"){
    regToExclude=paste(chr_regToInclude,":",mid_baitCoords-5000,"-",mid_baitCoords+5000,sep="")
    print(paste("regToExclude=",regToExclude))
}

chr_regToExclude=unlist(strsplit(regToExclude,split=":"))[1]
coord_regToExclude=as.numeric(unlist(strsplit(unlist(strsplit(regToExclude,split=":"))[2],"-")))
names(coord_regToExclude)=c("start","end")

data <- read.delim(fragsFile,skip=nskip,header=FALSE,sep=curSep)
colnames(data)=c("chr","start","end","score")
I.include=which(data[,"chr"]==chr_regToInclude & data[,"start"]<coord_regToInclude["end"] & data[,"end"]>coord_regToInclude["start"])
I.exclude=which(data[,"chr"]==chr_regToExclude & data[,"start"]<coord_regToExclude["end"] & data[,"end"]>coord_regToExclude["start"])
data[I.exclude,4]=NA

data.reg=data[I.include,]
mean.reg=mean(data.reg[,4],na.rm=TRUE)
print(mean.reg)


data.norm=data.frame(chr=data[,1], start=as.numeric(data[,2]), end=as.numeric(data[,3]), norm.score=as.numeric(data[,4]/mean.reg))
o <- order(data.norm[,1], data.norm[,2],data.norm[,3])
header=paste("track type=bedGraph name='",curName," normalised (factor=",mean.reg,")' description='",curName," normalised (factor=",mean.reg,")' visibility=full windowingFunction=maximum",sep="")
write.table(header,normFile,row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(data.norm[o,],file=normFile,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)

print("Done!")