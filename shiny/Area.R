cord = c("NDUFA13","NDUFB8","SDHA","UqCRC2","COX4+4L2","MTCO1","OSCP","VDAC1")#,"Dystrophin","DNA1")
mitochan = "VDAC1"
source("parseData.R", local = TRUE)

dat = getData("dat.txt",cord,mitochan)
dat = dat[dat$channel=="Area",]
boxplot(dat$value~dat$patient_id)
median(dat$value)

hipps = getData("C:\\GitHub\\mito_Bone\\dat.txt",cord,mitochan)
hipps = hipps[hipps$channel=="Area",]
boxplot(hipps$value~hipps$id)
median(hipps$value)