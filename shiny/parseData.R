library(data.table)

getData = function(fname,cord,mitochan="VDAC1"){
  dat = fread(fname,sep="\t",stringsAsFactors=FALSE,header=TRUE)
  dat$channel = gsub("GRIM19","NDUFA13",dat$channel)

  dat$ch = substring(dat$channel,regexpr("\\_[^\\_]*$", dat$channel)+1,nchar(dat$channel))
  dat = dat[dat$ch%in%cord,]

  dat$type = "Mean intensity"
  dat$type[grepl("LOG_",dat$channel)] = "Log mean intensity"
  dat$type[grepl("MED_",dat$channel)] = "Median intensity"
  dat$type[grepl("R_",dat$channel)] = "Ratio mean intensity (VDAC1)"
  dat$type[grepl("R_MED_",dat$channel)] = "Ratio median intensity (VDAC1)"
  dat$type[grepl("R_LOG_",dat$channel)] = "Ratio log mean intensity (VDAC1)"
  dat$type[grepl("Z_",dat$channel)] = "z-score"
  dat$outlier_diff = "NODIFF"
  dat$regression_diff = "NODIFF"
  dat$z_diff = "NODIFF"
  dat$z = 0

  dat$chstr = dat$ch
  transform = log
  dat_r = dat[dat$type=="Mean intensity",]
  dat_r$type = "r (VDAC1)"
  dat_theta = dat[dat$type=="Mean intensity",]
  dat_theta$type = "theta (VDAC1)"

  for(pid in unique(dat$patrep_id)){
   for(ch in unique(dat$ch)){
  	dt = dat[(dat$patrep_id==pid)&(dat$type=="Mean intensity"),]

	isch = as.character(dt$ch)==ch
	ismito = as.character(dt$ch)==mitochan
	prot = dt[isch,]
	mito = dt[ismito,]

	x = mito$value
	y = prot$value
	dat_r$value[(dat_r$patrep_id==pid)&(as.character(dat_r$ch)==ch)] = sqrt(x^2+y^2)
	dat_r$channel[(dat_r$patrep_id==pid)&(as.character(dat_r$ch)==ch)] = paste("RADIUS",ch,sep="_")
	dat_theta$value[(dat_theta$patrep_id==pid)&(as.character(dat_theta$ch)==ch)] = 360*atan(y/x)/(2*pi)
	dat_theta$channel[(dat_theta$patrep_id==pid)&(as.character(dat_theta$ch)==ch)] = paste("THETA",ch,sep="_")
    }
  }
  dat=rbind(dat,dat_r,dat_theta)
  return(dat)
}