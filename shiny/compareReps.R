cord = c("NDUFA13","NDUFB8","SDHA","UqCRC2","COX4+4L2","MTCO1","OSCP","VDAC1")#,"Dystrophin","DNA1")
mitochan = "VDAC1"
source("parseData.R", local = TRUE)

dat = getData("dat.txt",cord,mitochan)
patreps = unique(dat$patrep_id)

dats = lapply(patreps, function(x) updateDat(dat, "theta (VDAC1)", x, "NDUFB8", cord))
dat = do.call("rbind",dats)

# Uncommenting these lines would calculate correlations based on VDAC1 ratio instead
#dat = dat[!grepl("_LOG_",dat$channel),]
#dat = dat[grepl("R_",dat$channel),]

dat$ch = substring(dat$channel,regexpr("\\_[^\\_]*$", dat$channel)+1,nchar(dat$channel))
dat$chanid = paste(dat$patrep_id,dat$ch,sep="_")
agg = aggregate(dat,by=list(dat$chanid),FUN=mean)

r2 = agg$Group.1[agg$replicate==2]
r1 = gsub("R02","R01",r2)
r3 = gsub("R02","R03",r2)

agg = agg[agg$Group.1%in%c(r1,r2,r3),]
dat = dat[dat$chanid%in%c(r1,r2,r3),]

ids = agg$Group.1
ids = gsub("R01", "", ids)
ids = gsub("R02", "", ids)
ids = gsub("R03", "", ids)
agg$ids = ids

makePlot=function(xrep,yrep){
  xvals = agg$value[agg$replicate==xrep]
  yvals = agg$value[agg$replicate==yrep]
  ids = agg$ids[agg$replicate==xrep]
  cols = rep("darkgrey",length(ids))
  cols[substring(ids,1,1)=="P"] = "pink"
  axrng = range(c(xvals,yvals))
  axrng[2] = 1.1*axrng[2]
  corval = signif(cor(xvals,yvals),2)
  mlab = paste("Correlation:",corval)
  plot(xvals,yvals,type="n",xlab=paste("Mean theta, replicate",xrep),ylab=paste("Mean theta, replicate",yrep),xlim=axrng,ylim=axrng,main=mlab)
  abline(a=0,b=1,lwd=2,col="lightgrey")
  points(xvals,yvals,pch=16,col=cols)
  text(xvals,yvals,ids,pos = 4,cex=0.5)
}

pdf("CompareReps.pdf",width=10,height=10)
op=par(mfrow=c(2,2))
 makePlot(1,2)
 makePlot(1,3)
 makePlot(2,3)
par(op)
dev.off()

length(unique(dat$cell_id[dat$replicate==1]))
length(unique(dat$cell_id[dat$replicate==2]))
length(unique(dat$cell_id[dat$replicate==3]))

