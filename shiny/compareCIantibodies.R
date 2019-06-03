# Set working directory to byPatient
source("app.R")

unique(dat$type)
for(ty in c("Mean intensity","theta (VDAC1)","Ratio mean intensity (VDAC1)")){
  dt = dat[dat$type==ty]
  allres=list()
  for(pid in unique(dt$patient_id)){
    dp = dt[dt$patient_id==pid,]
    tres = t.test(dp$value[dp$ch=="NDUFB8"],dp$value[dp$ch=="GRIM19"])
    fres = list(tres$p.value,tres$conf.int[1],tres$conf.int[2])
    allres[[pid]]=fres
  }

  allres = data.frame(do.call(rbind,allres),stringsAsFactors=FALSE)
  colnames(allres)=c("p","lower_95","upper_95")
  allres$q = p.adjust(allres$p,"fdr")
  
  print(ty)
  print(allres)
}