source("dataFunctions.R")
source("parseData.R")

hypertest=function(x,m,y,n) {
    1-phyper(x-1,y,n-y,m)
  }

hiliteChannel = function(dat, hilite_ch = "NDUFB8", hilite_type = "theta (VDAC1)",alph=0.3){
 hcol = colorRamp(c("red","yellow","blue"),space="Lab")
 hcol_rgb = function(x, alpha=alph){
   vals = hcol(x)
   return(rgb(vals[1]/255,vals[2]/255,vals[3]/255,alpha))
 }

 dhilite = dat[(as.character(dat$ch)==hilite_ch)&(dat$type==ifelse(hilite_type=="2Dmito","theta (VDAC1)",hilite_type)),]
 dhilite$quant = ecdf(dhilite$value)(dhilite$value)
 
 dhcolours = sapply(dhilite$quant,hcol_rgb)
 names(dhcolours) = dhilite$cell_id
 return(dhcolours[dat$cell_id])
}

cord = c("NDUFB8","NDUFA13","SDHA","UqCRC2","COX4+4L2","MTCO1","OSCP","VDAC1","Dystrophin","DNA1")
cutcords = c(2.5,3.5,4.5,6.5,7.5,8.5)
cordlabs = c("CI","CII","CIII","CIV","CV","OMM","Cell")

complexes = c("CI","CI","CII","CIII","CIV","CIV","CV","Mito","Cell wall","Nucleus")
names(complexes) = cord

newids = c("P01","P02","P03","P04","P05","P06","P07","P08","P09","P10","C01","C02","C03")
names(newids) = c("P_M0838","P_M0517","P_M0966","P_M0284","P_M0694","P_m0917","P_M963","P_M0913","P_M1513","P_M0207","C_m1180","C_M1089","C_M1217")

revnew=names(newids)
names(revnew)=newids

dat = read.delim("mitocyto_merged_results with feature.csv",sep=",",stringsAsFactors=FALSE)

dat$Filename = gsub("M0207-10","M0207",dat$Filename)
dat$Filename = gsub("M0913-14","M0913",dat$Filename)
dat$Filename = gsub("M1180-16","m1180",dat$Filename)
dat$Filename = gsub("M1217-16","M1217",dat$Filename)
dat$Replicate = 1
colnames(dat) = c("value","id","channel","patient_type","patient_id","replicate")
dat$patient_type[dat$patient_id%in%c("m1180","M1217")]="control"
dat$patient_type[dat$patient_id%in%c("M0207","M0913")]="patient"

dat$patient_id = sub("-.*","",dat$patient_id)
#dat$cell_id = sub("-.*_","_",dat$cell_id)

# Specify which ids correspond to patients, and which to control
dat$patient_id = paste(toupper(substr(dat$patient_type,1,1)),dat$patient_id,sep="_")

dat$subject_group = "none"
sgroups = list(
Control = c("C_m1180","C_M1089","C_M1217"),
CI = c("P_M0838","P_M0517"), 
Deletion = c("P_M0966","P_M0284"),
"MT-TL1" = c("P_M0694","P_m0917","P_M963"),
"MT-TW" = c("P_M0207"),
"MT-TG" = c("P_M0913"),
"MT-TE" = c("P_M1513")
)
for(sgroup in names(sgroups)){
  dat$subject_group[dat$patient_id%in%sgroups[[sgroup]]] = sgroup
}

# More explicit group descriptions
subtext =c("Healthy control","Nuclear-encoded mutation\ncausing defect in complex I","Single, large-scale\nmtDNA deletion","Point mutation in mitochondrially\nencoded tRNA Leucine 1 (MT-TL1)","Point mutation in mitochondrially\nencoded tRNA (MT-TE)","Point mutation in mitochondrially\nencoded tRNA (MT-TG)","Point mutation in mitochondrially\nencoded tRNA (MT-TW)")
names(subtext) = c("Control", "CI", "Deletion", "MT-TL1", "MT-TE", "MT-TG", "MT-TW")

# Drop these subjects
#dat = dat[!dat$patient_id%in%c("P_m0164","P_M1144","P_M1600"),]
dat = dat[!dat$subject_group=="none",]

dat$patient_id = newids[dat$patient_id]
dat$patrep_id = paste(dat$patient_id,sprintf("R%02d",dat$replicate),sep="")

dat$cell_id = paste(dat$patrep_id,sprintf("%04d",dat$id),sep="_")
dat = dat[!duplicated(dat[,c("cell_id","channel")]),]

# Specify some colours for plotting
dat$colour = "black"
dat$colour[dat$patient_type=="patient"]=rgb(1,0,0,0.04)
dat$colour[dat$patient_type=="control"]=rgb(0,0,1,0.04)

# Drop these channels
todrop = c("103Rh","148Nd","151Eu","152S","154S","156Gd","165H","167Er","189Os","TOM22","DNA2","Dystrophin","DNA1")
cord = cord[!cord%in%todrop]
cutcords = cutcords[1:(length(cutcords)-1)]
cordlabs = cordlabs[1:(length(cordlabs)-1)]

channels = sort(unique(dat$channel))
channels = channels[sapply(channels,function(x) sum(sapply(todrop,function(y) grepl(y,x)))==0)]
dat = dat[dat$channel%in%channels,]

# Pretty names
getProtein=function(x){
 if(grepl("_",x)){
  substring(x,max(unlist(gregexpr(pattern="_",x)))+1,nchar(x))
 }else{
  x
 }
}
cns = dat$channel
cns = ifelse(grepl("DNA",cns),substring(cns,nchar(cns)-3, nchar(cns)),cns)
cns = substring(cns,unlist(gregexpr(pattern=" ",cns))+1,nchar(cns))
cns = sapply(cns, getProtein)
cns = ifelse(grepl("LOG_",dat$channel),paste("LOG",cns,sep="_"),cns)
cns = ifelse(grepl("MED_",dat$channel),paste("MED",cns,sep="_"),cns)#
cns = gsub("COX4","COX4+4L2",cns)
dat$channel = cns
channels = sort(unique(dat$channel))

dat = dat[!grepl(".png",dat$channel),]

dat$type = "mean intensity"
dat$type[grepl("LOG_",dat$channel)] = "log mean intensity"
dat$type[grepl("MED_",dat$channel)] = "median intensity"
dat$type[grepl("Area",dat$channel)] = "area"
dat$type[grepl("AspectRatio",dat$channel)] = "aspect ratio"
dat$type[grepl("Perimeter",dat$channel)] = "perimeter"
dat$type[grepl("Circularity",dat$channel)] = "circularity"

write.table(dat,"dat.txt",sep="\t",quote=FALSE,row.names=FALSE)

dat = as.data.frame(getData("dat.txt",cord,"VDAC1"))

#dat = updateDat(dat,"theta (VDAC1)",)

bychans = c("NDUFB8","SDHA","UqCRC2","COX4+4L2","OSCP","VDAC1")

cmax = max(dat$value[dat$channel%in%c("xCoord","yCoord")])

library(plotfunctions)
library(plotrix)
alph=1.0

pdf("SpatialReport.pdf",width=14,height=14)
for(bychan in bychans[bychans!="VDAC1"]){
  type = "theta (VDAC1)"
  if(bychan=="AspectRatio") type = "AspectRatio"
  if(bychan=="Area") type = "Area"
  if(bychan=="Perimeter") type = "Perimeter"
  if(bychan=="Circularity") type = "Circularity"
  bc = bychan
  if(type=="theta (VDAC1)") bc = paste("THETA",bychan,sep="_")

  dhcolours = hiliteChannel(dat,bychan,type,alp=1.0)
  hcol = colorRamp(c("red","yellow","blue"),space="Lab")
  hcol_rgb = function(x, alpha=alpha){
   vals = hcol(x)
   return(rgb(vals[1]/255,vals[2]/255,vals[3]/255,alpha))
  }

  for(pid in sort(unique(dat$patrep_id))){

    cexmax = 2
    cexmin = 0.3
    rmin = 5
    rmax = 50
    mitochan = "VDAC1"
    vals = dat$value[(as.character(dat$ch)==bychan)&(dat$type==ifelse(type=="2Dmito","theta (VDAC1)",type))]
    dat2 = as.data.frame(updateDat(dat,type,pid,bychan,cord),stringsAsFactors=FALSE)
    dat$outlier_diff[(dat$patrep_id==pid)&(dat$ch%in%unique(dat2$ch))]=as.character(dat2$outlier_diff)
    dat$regression_diff[(dat$patrep_id==pid)&(dat$ch%in%unique(dat2$ch))]=as.character(dat2$regression_diff)
    dat$z_diff[(dat$patrep_id==pid)&(dat$ch%in%unique(dat2$ch))]=as.character(dat2$z_diff)
    dt = dat[(as.character(dat$patrep_id)==pid),]
    
    cmax = max(dt$value[dt$channel%in%c("xCoord","yCoord")])

    N = length(unique(dt$id))

    difftype = "regression_diff"
    td = reshape(dt[,c("value","channel","cell_id")],idvar="cell_id",timevar="channel",direction="wide")
    tdiff = reshape(dt[,c(difftype,"channel","cell_id")],idvar="cell_id",timevar="channel",direction="wide")
    colnames(tdiff) = gsub(paste(difftype,".",sep=""),"",colnames(tdiff))
    colnames(td) = gsub("value.","",colnames(td))
    td$cex = cexmin+(td$Area-min(dat$value[dat$channel=="Area"]))/(max(dat$value[dat$channel=="Area"])-min(dat$value[dat$channel=="Area"]))*(cexmax-cexmin)
    td$rad = rmin+(td$Area-min(dat$value[dat$channel=="Area"]))/(max(dat$value[dat$channel=="Area"])-min(dat$value[dat$channel=="Area"]))*(rmax-rmin)

    bchan = bychan
    if(bychan%in%names(complexes)) bchan = paste(bchan,"(",complexes[bychan],")",sep="")
    mlab = paste(paste(pid,"coloured by",bchan,"N =",N),subtext[dt$subject_group[1]],sep="\n")

    # Deficiency gradient
    op=par(mfrow=c(2,2))
    plot(td$xCoord,max(td$yCoord)-td$yCoord,xlab="x-coordinate (px)",ylab="y-coordinate (px)",type="n",main=mlab,cex.lab=1.5,cex.axis=1.5,xlim=c(0,cmax),ylim=c(0,cmax))
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
    #points(td$xCoord,max(td$yCoord)-td$yCoord,pch=16,cex=td$cex,col=dhcolours[td$cell_id])
    symbols(td$xCoord,max(td$yCoord)-td$yCoord,td$rad,inches=FALSE,add=TRUE,fg=dhcolours[td$cell_id],bg=dhcolours[td$cell_id])
    gradientLegend(range(vals),color=sapply(ecdf(vals)(seq(min(vals),max(vals),length.out=10)),hcol_rgb),side=4,pos.num=4,pos=0.85,dec=2,n.seg=5)

    # Neighbourhood analysis
    dmat = as.matrix(dist(data.frame(x=td$xCoord,y=td$yCoord)))
    colnames(dmat)=td$cell_id
    rownames(dmat)=td$cell_id

    neighb = dmat<=5*mean(sqrt(td$Area/(2*pi)))
    sum(neighb)/length(td$cell_id)

    defect = tdiff[[bc]]=="BELOW"
    Ndef = sum(defect)
    Ntot = length(defect)

    pvals=rep(1,Ndef)
    samps = sample(seq_along(td$cell_id),4,replace=FALSE)

    neighbs = c()

    defects = td$cell_id[defect]
    for(j in seq_along(td$cell_id[defect])){
     i = td$cell_id[defect][j]
     nb = td$cell_id[neighb[i,]]
     nb = nb[nb!=i]
     neighbs = c(neighbs,nb)
     pvals[j] = hypertest(sum(neighb[i,]&defect),sum(neighb[i,]),Ndef-1,Ntot-1)
    }
    qvals = p.adjust(pvals,method="fdr")
    neighbs = unique(neighbs)
    defneighb = intersect(neighbs,defects)
    allhyp = hypertest(length(defneighb),length(neighbs),Ndef,Ntot)

    mlab2 = paste("No. defective:",Ndef,"No. in neighbourhood:",length(neighbs),"Overlap:",length(defneighb),"\np-value over-representation:",allhyp)
    plot(td$xCoord,max(td$yCoord)-td$yCoord,xlab="x-coordinate (px)",ylab="y-coordinate (px)",type="n",main=mlab2,cex.lab=1.5,cex.axis=1.5,xlim=c(0,cmax),ylim=c(0,cmax))
    symbols(td$xCoord,max(td$yCoord)-td$yCoord,td$rad,inches=FALSE,add=TRUE,fg="grey",bg="grey")

      #for(i in defects){
      # symbols(td$xCoord[neighb[i,]],max(td$yCoord)-td$yCoord[neighb[i,]],td$rad[neighb[i,]],inches=FALSE,add=TRUE,fg="blue",bg="blue")
      # symbols(td$xCoord[td$cell_id==i],max(td$yCoord)-td$yCoord[td$cell_id==i],td$rad[td$cell_id==i],inches=FALSE,add=TRUE,fg="red",bg="red")
      #}
       if(length(neighbs)>0) symbols(td$xCoord[td$cell_id%in%neighbs],max(td$yCoord)-td$yCoord[td$cell_id%in%neighbs],td$rad[td$cell_id%in%neighbs],inches=FALSE,add=TRUE,fg="blue",bg="blue")
       if(Ndef>0) symbols(td$xCoord[defect],max(td$yCoord)-td$yCoord[defect],td$rad[defect]/2,inches=FALSE,add=TRUE,fg="yellow",bg="yellow")


# RANDOM SAMPLING
    vals = dat$value[(as.character(dat$ch)==bychan)&(dat$type==ifelse(type=="2Dmito","theta (VDAC1)",type))]
    dat2 = as.data.frame(updateDat(dat,type,pid,bychan,cord),stringsAsFactors=FALSE)
    dat$outlier_diff[(dat$patrep_id==pid)&(dat$ch%in%unique(dat2$ch))]=as.character(dat2$outlier_diff)
    dat$regression_diff[(dat$patrep_id==pid)&(dat$ch%in%unique(dat2$ch))]=as.character(dat2$regression_diff)
    dat$z_diff[(dat$patrep_id==pid)&(dat$ch%in%unique(dat2$ch))]=as.character(dat2$z_diff)
    dt = dat[(as.character(dat$patrep_id)==pid),]

    uids = unique(dt$cell_id)
    mult = dim(dt)[1]/length(uids)
    sampind = sample(1:length(uids),length(uids),replace=FALSE)

    sampinds = rep(sampind,mult)
    sampinds = sampinds + rep((0:(mult-1))*length(sampind),each=length(sampind))
    dt$value = dt$value[sampinds]
    
        cmax = max(dt$value[dt$channel%in%c("xCoord","yCoord")])

    N = length(unique(dt$id))

    difftype = "regression_diff"
    td = reshape(dt[,c("value","channel","cell_id")],idvar="cell_id",timevar="channel",direction="wide")
    tdiff = reshape(dt[,c(difftype,"channel","cell_id")],idvar="cell_id",timevar="channel",direction="wide")
    colnames(tdiff) = gsub(paste(difftype,".",sep=""),"",colnames(tdiff))
    colnames(td) = gsub("value.","",colnames(td))
    td$cex = cexmin+(td$Area-min(dat$value[dat$channel=="Area"]))/(max(dat$value[dat$channel=="Area"])-min(dat$value[dat$channel=="Area"]))*(cexmax-cexmin)
    td$rad = rmin+(td$Area-min(dat$value[dat$channel=="Area"]))/(max(dat$value[dat$channel=="Area"])-min(dat$value[dat$channel=="Area"]))*(rmax-rmin)

    bchan = bychan
    if(bychan%in%names(complexes)) bchan = paste(bchan,"(",complexes[bychan],")",sep="")
    mlab = paste(paste(pid,"coloured by",bchan,"N =",N),subtext[dt$subject_group[1]],sep="\n")

    # Deficiency gradient
    plot(td$xCoord,max(td$yCoord)-td$yCoord,xlab="x-coordinate (px)",ylab="y-coordinate (px)",type="n",main=mlab,cex.lab=1.5,cex.axis=1.5,xlim=c(0,cmax),ylim=c(0,cmax))
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
    #points(td$xCoord,max(td$yCoord)-td$yCoord,pch=16,cex=td$cex,col=dhcolours[td$cell_id])
    symbols(td$xCoord,max(td$yCoord)-td$yCoord,td$rad,inches=FALSE,add=TRUE,fg=dhcolours[td$cell_id],bg=dhcolours[td$cell_id])
    gradientLegend(range(vals),color=sapply(ecdf(vals)(seq(min(vals),max(vals),length.out=10)),hcol_rgb),side=4,pos.num=4,pos=0.85,dec=2,n.seg=5)

    # Neighbourhood analysis
    dmat = as.matrix(dist(data.frame(x=td$xCoord,y=td$yCoord)))
    colnames(dmat)=td$cell_id
    rownames(dmat)=td$cell_id

    neighb = dmat<=5*mean(sqrt(td$Area/(2*pi)))
    sum(neighb)/length(td$cell_id)

    defect = tdiff[[bc]]=="BELOW"
    Ndef = sum(defect)
    Ntot = length(defect)

    pvals=rep(1,Ndef)
    samps = sample(seq_along(td$cell_id),4,replace=FALSE)

    neighbs = c()

    defects = td$cell_id[defect]
    for(j in seq_along(td$cell_id[defect])){
     i = td$cell_id[defect][j]
     nb = td$cell_id[neighb[i,]]
     nb = nb[nb!=i]
     neighbs = c(neighbs,nb)
     pvals[j] = hypertest(sum(neighb[i,]&defect),sum(neighb[i,]),Ndef-1,Ntot-1)
    }
    qvals = p.adjust(pvals,method="fdr")
    neighbs = unique(neighbs)
    defneighb = intersect(neighbs,defects)
    allhyp = hypertest(length(defneighb),length(neighbs),Ndef,Ntot)

    mlab2 = paste("No. defective:",Ndef,"No. in neighbourhood:",length(neighbs),"Overlap:",length(defneighb),"\np-value over-representation:",allhyp)
    plot(td$xCoord,max(td$yCoord)-td$yCoord,xlab="x-coordinate (px)",ylab="y-coordinate (px)",type="n",main=mlab2,cex.lab=1.5,cex.axis=1.5,xlim=c(0,cmax),ylim=c(0,cmax))
    symbols(td$xCoord,max(td$yCoord)-td$yCoord,td$rad,inches=FALSE,add=TRUE,fg="grey",bg="grey")

      #for(i in defects){
      # symbols(td$xCoord[neighb[i,]],max(td$yCoord)-td$yCoord[neighb[i,]],td$rad[neighb[i,]],inches=FALSE,add=TRUE,fg="blue",bg="blue")
      # symbols(td$xCoord[td$cell_id==i],max(td$yCoord)-td$yCoord[td$cell_id==i],td$rad[td$cell_id==i],inches=FALSE,add=TRUE,fg="red",bg="red")
      #}
       if(length(neighbs)>0) symbols(td$xCoord[td$cell_id%in%neighbs],max(td$yCoord)-td$yCoord[td$cell_id%in%neighbs],td$rad[td$cell_id%in%neighbs],inches=FALSE,add=TRUE,fg="blue",bg="blue")
       if(Ndef>0) symbols(td$xCoord[defect],max(td$yCoord)-td$yCoord[defect],td$rad[defect]/2,inches=FALSE,add=TRUE,fg="yellow",bg="yellow")


    par(op)
  }
}
dev.off()


