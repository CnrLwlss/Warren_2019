dat=read.delim("Warren_2020.csv",sep=",",stringsAsFactors=FALSE)
getprot = function(ch) ifelse(grepl(" ",ch),strsplit(ch," ")[[1]][2],ch)
dat$Channel=sapply(dat$Channel,getprot)
write.table(dat,"Warren_2020_chans.csv",sep=",",quote=FALSE,row.names=FALSE)