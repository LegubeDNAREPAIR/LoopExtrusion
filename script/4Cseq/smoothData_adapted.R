#smoothData_adapted.R : adapted from smoothData (https://github.com/bbcf/bbcfutils/blob/master/R/smoothData.R) to be used in a snakemake pipeline

.libPaths(c(.libPaths(),"/home/rochevin/R-packages/"))
setwd(snakemake@params[["dir"]]) ;
require(dplyr)
options(stringsAsFactors = FALSE)
options(scipen=999) # disable the scientific notation


fragsFile <- snakemake@input[[1]]
nFragsPerWin <- snakemake@params[["nFragsPerWin"]] %>% as.numeric()
curName <- snakemake@params[["curName"]]
resFile <- snakemake@output[[1]]
regToExclude <- snakemake@params[["regToExclude"]]
scoreCol <- snakemake@params[["scoreCol"]] %>% as.numeric()

skipN=1

print(paste("fragsFile=",fragsFile))
print(paste("nFragsPerWin=",nFragsPerWin))
print(paste("curName=",curName))
print(paste("resFile=",resFile))
print(paste("regToExclude=",regToExclude))
print(paste("ScoreColumn from input=",scoreCol))


updateFactors <- function(myData)
{
    myData[] <- lapply(myData, function(x) x[,drop=TRUE])
    return(myData)
}

#getRegions <- function(myData,baitCoord,interRegCoord)
getRegions <- function(myData,interRegCoord)
{
    i_baitChr <- which(myData[,1] == interRegCoord[1])
    notBaitChr <- updateFactors(myData[which(myData[,1] != interRegCoord[1]),])
    ##i <- which(myData[i_baitChr,3] < as.numeric(interRegCoord[2]))
    i <- which(myData[i_baitChr,3] < as.numeric(interRegCoord[2])) ##fragments at the border are excluded
    regUp <- updateFactors(myData[i_baitChr[i],])
    #i <- which(myData[i_baitChr,3] > as.numeric(interRegCoord[3]))
    i <- which(myData[i_baitChr,2] > as.numeric(interRegCoord[3]))
    regDown <- updateFactors(myData[i_baitChr[i],])
    return(list(othersChrs=notBaitChr,upstream_of_interactiveReg=regUp,down_of_interactiveReg=regDown))
}

mergeLists <- function(x,y)
{
    if (length(x) == 0){return(y)}
    if (length(y) == 0){return(x)}
    i = intersect(names(y), names(x))
    j <- setdiff(names(x),names(y))
    k <- setdiff(names(y),names(x))
    z <- list()
    if(length(i)>0){for(n in i){z[[n]] <- rbind(x[[n]],y[[n]])}}
    if(length(j)>0){for(n in j){z[[n]]<-x[[n]]}}
    if(length(k)>0){for(n in k){z[[n]]<-y[[n]]}}
    return(z)
}


splitData <- function(myData)
{
    z <- list();
    for(i in names(table(myData[,1]))){j <- which(myData[,1]==i);z[[i]]<-cbind(myData[j,],rownames(myData)[j])}
    return(z)
}


groupFrag_betweenIndices <- function(myVector,s,e)
{
    return(lapply(1:length(s),function(j){mean(myVector[s[j]:e[j]]) }))
}


generateSmoothedWindows <- function(myData,n)
{
    data_byChr <- splitData(myData)
    if(length(data_byChr)>0){
        newWindows <- lapply(1:length(data_byChr),
                             function(i){if(nrow(data_byChr[[i]])>0){
                                 s=c();e=c();
                                 mid=round(n/2);
                                 nprev=mid;nnext=mid;
                                 #if((n/2)==round(n/2)){nprev=mid-1;nnext=mid}else{nprev=mid;nnext=mid}
                                 mid=seq(1,nrow(data_byChr[[i]]),by=1)
                                 istart=mid-rep(nprev,length(mid))
                                 iend=mid+rep(nnext,length(mid))
                                 ipairs <- cbind(istart,iend,mid)
                                 ipairs[which(ipairs[,2]>length(mid)),2]=length(mid)
                                 ipairs[which(ipairs[,1]<0),1]=1
                                 
                                 cbind(data_byChr[[i]][mid,2],data_byChr[[i]][mid,3],groupFrag_betweenIndices(data_byChr[[i]][,4],ipairs[,1],ipairs[,2]))
                             } #end if nrow>0
                             } #end function(i)
                             
        )
        names(newWindows)=names(data_byChr)}
    else{newWindows=c()}
    return(newWindows)
}

stringsAsFactors=FALSE
print(paste("read file ",fragsFile,sep=""))
allFrags <- read.delim(fragsFile,header=FALSE,skip=skipN)
data <- data.frame("chr"=as.factor(allFrags[,1]),"start"=as.integer(allFrags[,2]),"end"=as.integer(allFrags[,3]),"score"=as.numeric(allFrags[,scoreCol]),stringsAsFactors=FALSE)
print("dim(data)")
print(dim(data))
if(nrow(data)==0){print("Warning!! file is empty!!")}
allSmoothedWindows <- generateSmoothedWindows(data,nFragsPerWin)

if(nchar(regToExclude)>2)
{
    regToExcludeSplit <- unlist(strsplit(gsub("-",":",regToExclude,perl=TRUE),":"))
    data_splitted <- getRegions(data,regToExcludeSplit)
    othersChrs_allWindows <- generateSmoothedWindows(data_splitted[[1]],nFragsPerWin)
    upInteractiveReg_allWindows <- generateSmoothedWindows(data_splitted[[2]],nFragsPerWin)
    downInteractiveReg_allWindows <- generateSmoothedWindows(data_splitted[[3]],nFragsPerWin)
    allSmoothedWindows <- mergeLists(othersChrs_allWindows,upInteractiveReg_allWindows)
    allSmoothedWindows <- mergeLists(allSmoothedWindows,downInteractiveReg_allWindows)
    header=paste("track type=bedGraph name='",curName," (all smoothed windows - ",nFragsPerWin," fragments per window)' description='",curName," (all smoothed windows - ",nFragsPerWin," fragments per window - region excluded=",regToExclude,")' visibility=full windowingFunction=maximum",sep="")
}

if(nchar(regToExclude)<=2)
{
    header=paste("track type=bedGraph name='",curName," (all smoothed windows - ",nFragsPerWin," fragments per window)' description='",curName," (all smoothed windows - ",nFragsPerWin," fragments per window)' visibility=full windowingFunction=maximum",sep="")
}



print("Will write file")
write.table(header,resFile,row.names=FALSE,col.names=FALSE,quote=FALSE)
#lapply(1:length(allSmoothedWindows),function(i){write.table(cbind(rep(names(allSmoothedWindows)[i],nrow(allSmoothedWindows[[i]])),allSmoothedWindows[[i]]),resFile,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE);return(nrow(allSmoothedWindows[[i]]))})
outlapply=lapply(1:length(allSmoothedWindows),function(i){j=which(allSmoothedWindows[[i]][,3]>0);write.table(cbind(rep(names(allSmoothedWindows)[i],length(j)),allSmoothedWindows[[i]][j,]),resFile,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE);return(nrow(allSmoothedWindows[[i]][j,]))})

print("*****************")
print(paste("Smoothing of ",fragsFile," done!"))
print(paste("resfile=",resFile))
print("*****************")
