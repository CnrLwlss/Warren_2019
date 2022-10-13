mitoplot = function(dat,pat,chan,pisd=qnorm(0.975),lab_inner="",
                    summary="Mean intensity",xlim=NA,ylim=NA,
                    makeplot=TRUE,makelinreg=TRUE,makectrl=TRUE,makepat=TRUE,makeclass=TRUE,
                    bootstrap=FALSE,seed = FALSE,mitochan="VDAC1",patalpha = 0.25,ctrlalpha=0.15,hilite_ids=c(),clust_class=c()){
                    
 colnames(dat) = tolower(colnames(dat))

 xlab = paste0("log(",mitochan,")")
 ylab = paste0("log(",chan,")")

 xpat = log(dat$value[(dat$ch==mitochan)&(dat$fn==pat)&(dat$type==summary)])
 ypat = log(dat$value[(dat$ch==chan)&(dat$fn==pat)&(dat$type==summary)])
 keepind = (!is.na(xpat))&(!is.na(ypat))
 xpat = xpat[keepind]
 ypat = ypat[keepind]
 yids = dat$id[(dat$ch==chan)&(dat$fn==pat)&(dat$type==summary)][keepind]
 xctrl = log(dat$value[(dat$ch==mitochan)&(grepl("C",dat$fn))&(dat$type==summary)])
 yctrl = log(dat$value[(dat$ch==chan)&(grepl("C",dat$fn))&(dat$type==summary)])
 xctrl = xctrl[(!is.na(xctrl))&(!is.na(yctrl))]
 yctrl = yctrl[(!is.na(xctrl))&(!is.na(yctrl))]

 if(bootstrap){
   if(typeof(seed)!="logical") set.seed(seed)
   patindex = sample(1:length(xpat),replace=TRUE)
   ctrlindex = sample(1:length(xctrl),replace=TRUE)
   xpat = xpat[patindex]
   ypat = ypat[patindex]
   xctrl = xctrl[ctrlindex]
   yctrl = yctrl[ctrlindex]
 }

 if(length(xlim)!=2) xlim = range(c(xpat,xctrl),finite=TRUE,na.rm=TRUE)
 if(length(ylim)!=2) ylim = range(c(ypat,yctrl),finite=TRUE,na.rm=TRUE)
 xsyn = seq(xlim[1],xlim[2],length.out=500)
 
 N = length(xpat)
 
 mod = lm(yctrl~xctrl)
 
 # https://en.wikipedia.org/wiki/Prediction_interval 
 pred_syn = predict(mod, newdata = data.frame(xctrl=xsyn), se.fit=TRUE, interval = "prediction",na.action=na.omit, level=0.95)$fit
 mid_syn = pred_syn[,'fit']
 up_syn = pred_syn[,'upr']
 low_syn = pred_syn[,'lwr']
 psd_syn = pisd*(up_syn - low_syn)/(2*qnorm(0.975))
 up_syn = mid_syn + psd_syn
 low_syn = mid_syn - psd_syn
 
 pred_dat = predict(mod, newdata=data.frame(xctrl=xpat), se.fit=TRUE,  interval = "prediction",na.action=na.omit, level=0.95)$fit
 mid_dat = pred_dat[,'fit']
 up_dat = pred_dat[,'upr']
 low_dat = pred_dat[,'lwr']
 psd_dat = pisd*(up_dat - low_dat)/(2*qnorm(0.975))
 up_dat = mid_dat + psd_dat
 low_dat = mid_dat - psd_dat

 z = (ypat - mid_dat)/((up_dat - low_dat)/(2*qnorm(0.975)))

 N = length(xpat)
 class = rep("Medium",N)
 class[ypat>up_dat] = "High"
 class[ypat<low_dat] = "Low"

 if(length(clust_class)!=0){
  if(length(clust_class)==N){
   class=clust_class
  }else{
   print("Incorrect number of classes in clust_class!")
   return()
  }
 }

 res =list(med = sum(class=="Medium"), low = sum(class=="Low"), high = sum(class=="High"), N = N, z = z,class=class)
 
 if(makeplot){
  ccol = rep(rgb(0,0,1,patalpha),N)
  ccol[class=="High"]=rgb(0,1,0,patalpha)
  ccol[class=="Low"]=rgb(1,0,0,patalpha)
  
  mlab = paste0(pat,"\n")
  if(makelinreg) mlab = paste0(pat,"\n Pred. int: ",signif(pisd,3))
  if(makepat) mlab = paste0(pat,"\nN: ",N)
  if(makeclass) mlab = paste0(pat,"\nN: ",N," Deficient: ",signif(100*res$low/res$N,3),"%")

  plot(NULL,pch=16,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,cex.axis=1.55,cex.lab=2.55,main=mlab,cex.main=1.9)
  if(makelinreg){
   points(xsyn,up_syn,lwd=4,lty=2,type="l")
   points(xsyn,mid_syn,lwd=6,lty=1,type="l")
   points(xsyn,low_syn,lwd=4,lty=2,type="l")
  }
  if(makectrl) points(xctrl,yctrl,col=rgb(0,0,0,ctrlalpha),cex=2,pch=16)
  if(makeclass){ccol2=ccol}else{ccol2 = rgb(1,0,1,patalpha)}
  if(makepat) {
   points(xpat,ypat,col=ccol2,cex=2,pch=16)
   points(xpat[yids%in%hilite_ids],ypat[yids%in%hilite_ids],cex=2.1,col=rgb(0,0,0),lwd=2)
  }
  text(x=xlim[1],y=ylim[2],lab_inner,cex=2,adj=c(0, NA))
 }
 return(res)
}

labelpanel = function(txt,cex=2.5){
 # https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/
 op2 = par(xpd=NA)
 di = dev.size("in")
 x <- grconvertX(c(0, di[1]), from="in", to="user")
 y <- grconvertY(c(0, di[2]), from="in", to="user")
 fig <- par("fig")
 x <- x[1] + (x[2] - x[1]) * fig[1:2]
 y <- y[1] + (y[2] - y[1]) * fig[3:4]
 x <- x[1] + strwidth(txt, cex=cex) / 2
 y <- y[2] - strheight(txt, cex=cex) / 2
 text(x, 0.95*y, txt, cex=cex,adj=c(0.5, NA))
 par(op2)
}

# COMBS
Nchans = 3
