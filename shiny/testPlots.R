source("plotFunctions.R", local = TRUE)
source("dataFunctions.R", local = TRUE)

subtext =c("healthy control","nuclear-encoded mutation in CI","single, large-scale mtDNA deletion","point mutation in mito. encoded tRNA Leucine 1 (MT-TL1)","point mutation in mito. encoded tRNA (MT-TE)","point mutation in mito. encoded tRNA (MT-TG)","point mutation in mito. encoded tRNA (MT-TW)")
names(subtext) = c("Control", "CI", "Deletion", "MT-TL1", "MT-TE", "MT-TG", "MT-TW")

cutcords = c(2.5,3.5,4.5,6.5,7.5)#,8.5)
cordlabs = c("CI","CII","CIII","CIV","CV","OMM")#,"Cell")
cord = c("NDUFA13","NDUFB8","SDHA","UqCRC2","COX4+4L2","MTCO1","OSCP","VDAC1")#,"Dystrophin","DNA1")
chlabs = c("CI","CI","CII","CIII","CIV","CIV","CV","OMM")
names(chlabs) = cord
mitochan = "VDAC1"

source("parseData.R", local = TRUE)

fulldat = "dat.txt"

alldat = getData(fulldat,cord,mitochan)

dat = updateDat(alldat[(alldat$patrep_id=="P05R01")&(alldat$type=="Mean intensity"),],"Mean intensity","P05R01","NDUFB8",cord)
cdat = alldat[(alldat$subject_group=="Control")&(alldat$type=="Mean intensity"),]
rdat =  alldat[(alldat$patrep_id=="P05R01")&(alldat$type=="theta (VDAC1)"),]
rdat$hcol = hiliteChannel(rdat, "NDUFB8", "theta (VDAC1)",alph=0.3)

dat = dat[!duplicated(dat[,c("cell_id","channel","type")]),]

xvals = dat$value[dat$channel=="VDAC1"]
yvals = dat$value[dat$channel=="NDUFB8"]

rawlim = c(0,max(c(xvals,yvals)))
loglim = c(0,max(log(c(xvals,yvals))))

pdf("AmyTest.pdf",width=12,height=6)
op=par(mfrow=c(1,2))
#plot(xvals,yvals,xlim=rawlim,ylim=rawlim)
plot(log(xvals),log(yvals),xlim=loglim,ylim=loglim,pch=16,col=rgb(0,0,0,0.3))

dat$cluster=1
arrayplot(dat,cdat,rdat,cord,"NDUFB8",hichan="NDUFB8")

par(op)
dev.off()
