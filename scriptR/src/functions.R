#Functions

loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}
get_sub_matrix_Hic <- function(file,bed,binsize){
    HTC <- loadRData(file)
    rows.chr <- x_intervals(HTC)
    val1 <- HTC %>% intdata %>% sum
    # HTC <- normPerReads(HTC)
    
    overlaps <- as(findOverlaps(bed,rows.chr), "List")
    overlaps.bin <- lapply(overlaps,length) %>% unlist()
    overlaps <- overlaps[overlaps.bin == binsize]
    
    res <- lapply(1:length(overlaps),function(i){
        x <- overlaps[[i]]
        if(as.character(strand(bed[i])) == "-"){
            x <- rev(x)
        }
        intdata(HTC)[x,x]
    })
    return(list(res,val1))
}
get_nb_read_by_chr <- function(file){
    HTC <- loadRData(file)
    HTC %>% intdata %>% sum
}