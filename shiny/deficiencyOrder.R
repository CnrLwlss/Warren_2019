install.packages("klaR")
library(klaR)
library(plyr)
# Set working directory to byPatient
source("app.R")

countFibres = function(dat, pat, cord, type="theta (VDAC1)", difftype="regression_diff", position = "BELOW", clevel = 0.95){
  dt = data.frame(updateDat(dat, type , pat, "NDUFB8", cord, clevel))
  Nfibres = length(unique(dt$cell_id[dt$patrep_id==pat]))
  formul = paste(difftype,"ch",sep="~")
  res = data.frame(aggregate(eval(parse(text=formul)),data=dt,function(x) sum(x==position)),stringsAsFactors=FALSE)
  colnames(res) = c("ch","count")
  res = rbind(res, data.frame(ch="TOTAL",count=Nfibres))
  return(res)
}

cats = function(pat,dat,type,difftype,chans, clevel = 0.95){
 dt = data.frame(updateDat(dat, type , pat, "NDUFB8", chans, clevel))
 Nfibres = length(unique(dt$cell_id[dt$patrep_id==pat]))
 wtab = reshape(dt[,c(difftype,"ch","cell_id")],idvar="cell_id",timevar="ch",direction="wide")
 colnames(wtab)=gsub(paste(difftype,".",sep=""),"",colnames(wtab))
 wtab = wtab[,c("cell_id",chans)]
 return(wtab)
}

countCategories = function(pat,dat,type,difftype,chans,clevel=0.95){
 wtab = cats(pat,dat,type,difftype,chans,clevel)
 combs = data.frame(unique(wtab[,2:length(colnames(wtab))]))
 rownames(combs) = 1:length(combs[,1])
 corder = do.call(order, as.list(wtab))

 wdat = wtab[,2:length(colnames(wtab))]
 colnames(wdat)= substr(colnames(wdat),1,4)
 wcount = count(wdat, vars=names(wdat))
 colnames(wcount) = c(chans,"freq")
 wcount$percent = signif(100*wcount$freq/length(wtab$cell_id),2)
 return(wcount[order(wcount$freq,decreasing=TRUE),])
}

chs = unique(dat$ch)
pats = unique(dat$patrep_id)

getcount = function(x) as.numeric(countFibres(dat,x,cord,type="theta (VDAC1)", difftype="regression_diff", position = "BELOW",clevel=0.99)$count)
ctab = sapply(pats,getcount)
rownames(ctab) = c(cord,"TOTAL")
ctab = data.frame(t(ctab))
ctab$VDAC1=NULL
ptab = signif(100*sweep(ctab,1,ctab$TOTAL,"/"),2)
ptab$TOTAL=NULL

neworder = sort(colSums(ptab),decreasing=TRUE)
ptab = ptab[,names(neworder)]
ctab = ctab[,c(names(neworder),"TOTAL")]

chans = names(neworder)
chans = c("NDUFB8", "MTCO1", "UqCRC2")
#chans = c("NDUFB8", "MTCO1", "UqCRC2","SDHA","OSCP")
chans = gsub("\\.","+",chans)

for(pat in unique(dat$patrep_id)){
 print(pat)
 print(countCategories(pat,dat,type="theta (VDAC1)", difftype="regression_diff",chans))
}

wtab = cats("P03R01",dat,type,difftype,chans)
wtab = wtab[,2:length(colnames(wtab))]
kmodes(wtab, 5, iter.max = 1000, weighted = FALSE )



