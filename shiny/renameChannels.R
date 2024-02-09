dat=read.delim("Warren_2020.csv",sep=",",stringsAsFactors=FALSE)
metals = c("153Eu_153Eu ","158Gd_158Gd ","160Gd_160Gd ","161Dy_161Dy ","164Dy_164Dy ","166Er_166Er ","168Er_168Er ","172Yb_172Yb ","174Yb_174Yb ","176Yb_176Yb ")

getprot = function(ch,metals) {
  for(met in metals) ch = gsub(met,"",ch)
  return(ch)
}

dat$Channel=sapply(dat$Channel,getprot,metals)
write.table(dat,"Warren_2020_chans.csv",sep=",",quote=FALSE,row.names=FALSE)