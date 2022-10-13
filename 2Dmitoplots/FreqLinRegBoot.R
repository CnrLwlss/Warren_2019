source("calcPlotDefect.R")
dat = read.delim("../shiny/dat.txt",stringsAsFactors=FALSE,sep="\t")
dat$ch = dat$channel
dat$fn = dat$patient_id

cord = c("GRIM19","NDUFB8","SDHA","UqCRC2","COX4+4L2","MTCO1","OSCP","VDAC1")#,"Dystrophin","DNA1")
chlabs = c("CI","CI","CII","CIII","CIV","CIV","CV","OMM")
names(chlabs) = cord
mitochan = "VDAC1"

subtext =c("healthy control","nuclear-encoded mutation in CI","single, large-scale mtDNA deletion","point mutation in mito. encoded tRNA Leucine 1 (MT-TL1)","point mutation in mito. encoded tRNA (MT-TE)","point mutation in mito. encoded tRNA (MT-TG)","point mutation in mito. encoded tRNA (MT-TW)")
names(subtext) = c("Control", "CI", "Deletion", "MT-TL1", "MT-TE", "MT-TG", "MT-TW")

subjs = sort(unique(dat$patient_id))
ctrls = subjs[grep("C",subjs)]
pats = subjs[grep("P",subjs)]

dir.create("BootstrapParticles", showWarnings = FALSE)
dir.create("2Dmitoplots", showWarnings = FALSE)

getprops = function(dat,pat,chan,summary="mean intensity",mitochan="VDAC1"){
  bs = mitoplot(dat,pat,chan,summary=summary,mitochan=mitochan,bootstrap=TRUE,makeplot=FALSE)
  props = c(bs$low/bs$N,bs$med/bs$N,bs$high/bs$N)
  names(props) = c("low","med","high")
  return(props)
}

for(pat in pats){
 for(chan in cord[cord!=mitochan]){
  print(paste(pat,chan))
  png(file.path("2Dmitoplots",paste0(pat,"_",chan,".png")),width=1500,height=1500,pointsize=30)
  op = par(mar= c(4.4, 5.2, 4, 0.1))
  mitoplot(dat,pat,chan,summary="mean intensity",mitochan=mitochan,lab_inner=chlabs[chan])
  par(op)
  props = t(replicate(1000,getprops(dat,pat,chan)))
  write.table(props,file=file.path("BootstrapParticles",paste0(pat,"_",chan,".txt")),sep="\t",quote=FALSE,row.names=FALSE)  
  dev.off()
 }
}
