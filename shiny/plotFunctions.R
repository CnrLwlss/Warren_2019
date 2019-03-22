getwide = function(dat, column = "value"){
 dw = reshape(subset(dat,select=c("cell_id","ch",column)),idvar="cell_id",timevar="ch",direction="wide")
 colnames(dw) = gsub(paste(column,".",sep=""),"",colnames(dw))
 return(dw)
}

stripAlpha = function(x) do.call(rgb,as.list(as.numeric(col2rgb(x)/255))) # Strip alpha channel

schart = function(dsig,csig,ids,subtext,nums,cord,cutcords,cordlabs,subjectLabel="",hichan=" ",axrngCheck = TRUE, axrng=c(-5,5),jitpat = "jitr",jitctrl="jitl",showControls=TRUE){
	dtype = unique(dsig$type)[1]
	
    #if(length(ids)>nmax) {ids = sample(ids,nmax)}
	hl = dsig[dsig$cell_id%in%ids,]
    hl = hl[order(hl$jit),]
    subjgrp = unique(dsig$subject_group)
    glab=paste("Group: ",subjgrp,", Subject: ",subjectLabel,sep="")
	description = subtext[subjgrp]
	substring(description,1,1)=toupper(substring(description,1,1))
	glab = paste(description,subjectLabel,sep="\n")
    main=paste(glab,", N = ",length(unique(dsig$id)),", Selected = ",length(unique(hl$cell_id)),sep="")

    numclust = max(dsig$cluster)
    rcols = rainbow(numclust, alpha = 0.2)
    if(hichan == " "){
      cols = rcols[dsig$cluster]
      bcol = rgb(0.89,0.89,0.89)
      hicol = rgb(0,0,0,0.75)
    }else{
      #cols = rgb(0,0,0,0.05)
	  cols = dsig$hcol
      bcol = rgb(1,1,1)
	  #hicol = rgb(1,0,0,0.3)
	  hicol = rgb(0,0,0,0.75)
    }

    op = par(mar=c(6, 4, 2, 0) + 0.1)
    #if(grepl("Ratio",dtype)){logval="y"}else{logval=""}
	if(grepl("z-score",dtype)){logval=""}else{logval="y"}
	if(grepl("theta",dtype)) logval=""
	if(axrngCheck){ 
	  axrng = range(c(dsig$value,csig$value))
	  if(grepl("theta",dtype)) axrng = c(0,90)
	}
    plot(NULL,xlab="",ylab = dtype,main = main,cex.main=0.75,las=2,axes=FALSE,log=logval,type="n",ylim=axrng,xlim=c(min(dsig$jitl),max(dsig$jitr)))
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bcol)  

	grid(NA,NULL,equilogs=FALSE,col="black",lwd=1.5)
	#if(showCentres) {for(i in 1:numclust) {
	#  yvals = vals$km$centers[i,sapply(cord,agrep,colnames(vals$km$centers))]
	#  xvals = 1:length(yvals) - 0.25*jwidth
	#  points(xvals, yvals, type="l", col=rcols[i], lwd=3)
	#  }  
	#}
    dsig$cex=0.3
	dsig$cex[dsig$outlier_diff%in%c("ABOVE","BELOW")]=0.5
    points(dsig[[jitpat]],dsig$value,col = cols, pch = 16, cex=dsig$cex)
	if(showControls) points(csig[[jitctrl]],csig$value,col=rgb(0,0,0,0.1),pch=16,cex=0.3)
    axis(2)
    axis(1,at=1:length(nums),labels=names(nums),las=2)
    abline(v = cutcords,lwd=1.5,lty=2,col="darkgrey")
    text(rowMeans(embed(c(0.5,cutcords,length(cord)+0.5),2)),max(axrng),cordlabs,pos=1,cex=0.75)
    for(i in ids){
	  newcols = sapply(hl$hcol[hl$cell_id==i],stripAlpha)
	  points(hl[[jitpat]][hl$cell_id==i],hl$value[hl$cell_id==i],type="b",lwd=0.5,col=newcols)
    }
    par(op)
}  
	
arrayplot = function(pimc,cimc,rdat,cord,ch,ids=c(),reg_diff=c(),mitochan="VDAC1",hichan=" ",showControls=TRUE,axrngCheck = TRUE, axrng=c(-5,5),chlab="", logify = TRUE){
	 pid = unique(pimc$patient_id)
	 ctrlcol = rgb(0,0,0,0.15)
	 patcol = rgb(1,0,0,0.25)
	 if(logify){transform = log}else{transform = identity}
	 patvals = transform(pimc$value)
     ctrlvals = transform(cimc$value)
	 if(axrngCheck) axrng = range(c(patvals,ctrlvals),na.rm=TRUE)
	 
	numclust = max(pimc$cluster)
    rcols = rainbow(numclust, alpha = 0.2)
	
	cols = "red"
	colsel = "red"
    if(hichan == " "){
      cols = rcols[pimc$cluster]
      bcol = rgb(0.89,0.89,0.89)
      hicol = rgb(0,0,0,0.75)
    }else{
      #cols = rgb(0,0,0,0.05)
	  ratd = rdat[as.character(rdat$ch)==as.character(ch)]
	  cols = ratd$hcol
	  colsel = sapply(ratd$hcol[ratd$cell_id%in%ids], stripAlpha)
      bcol = rgb(1,1,1)
	  #hicol = rgb(1,0,0,0.3)
	  hicol = rgb(0,0,0,0.75)
    }
	 
     #for(ch in cord[cord!=mitochan]){
	  regression_diff = reg_diff[paste(ch,pimc$cell_id[as.character(pimc$ch)==ch])]

      xctrl = transform(cimc$value[(as.character(cimc$ch)==mitochan)&(cimc$subject_group=="Control")])
      yctrl = transform(cimc$value[(as.character(cimc$ch)==ch)&(cimc$subject_group=="Control")])
	  
	  pimcx = pimc[as.character(pimc$ch)==mitochan,]
	  getch = as.character(pimc$ch)==ch # Why the %$£& do I have to do this to make the filter work here?!
	  pimcy = pimc[getch,]
	  
	  xvals = transform(pimcx$value)
      yvals = transform(pimcy$value)
	  
	  xsel = xvals[pimcx$cell_id%in%ids]
	  ysel = yvals[pimcy$cell_id%in%ids]


	  N = length(xvals)

      xsyn = seq(min(axrng),max(axrng),length.out=50)
      mod = lm(yctrl~xctrl)
      pred = predict(mod,newdata = data.frame(xctrl=xsyn), se.fit=TRUE,  interval = "prediction",na.action=na.omit)$fit
	  mid = pred[,1]
      up = pred[,3]
      low = pred[,2]
      psd = (up - low)/(2*1.96)
      upz = mid+3*psd
      lowz = mid-3*psd
      #upy = ifelse(length(unique(up))>1,approxfun(xsyn,up),function(x) 0)
      #lowy = ifelse(length(unique(up))>1,approxfun(xsyn,low),function(x) 0)
	  upy = approxfun(xsyn,up)
      lowy = approxfun(xsyn,low)
	  
	  normcex = 0.75
	  cexvec = rep(normcex,N)
	  cexvec[regression_diff=="BELOW"|regression_diff=="ABOVE"]=1.125
	  cexvecsel = 1.5*cexvec[pimcx$cell_id%in%ids]
      frac = (sum(regression_diff=="BELOW")+sum(regression_diff=="ABOVE"))/N
      mlab = paste("N =",N,"Above =",sum(regression_diff=="ABOVE"),"Below =",sum(regression_diff=="BELOW"))#,"F =",signif(frac,2))
      #mlab = paste(paste(pid,unique(pimc$subject_group)),mlab,sep="\n")
	  #whatlog = ifelse(identical(log,transform),"","xy")
	  whatlog = ""
	  if(logify){
	    xlab=paste("log(",mitochan,")",sep="")
	    ylab=paste("log(",ch,")",sep="")
	  }else{
	    xlab = mitochan
		ylab = ch
	  }
      plot(NULL,xlab=xlab,ylab=ylab,main=mlab,type="n",xlim=axrng,ylim=axrng,log=whatlog,cex.main=1.0,cex.lab=1.5)
      points(xsyn,mid,type="l",lty=1,lwd=2,col=ctrlcol)
      points(xsyn,up,type="l",lty=2,lwd=2,col=ctrlcol)
      points(xsyn,low,type="l",lty=2,lwd=2,col=ctrlcol)
      #points(xsyn,lowz,type="l",lty=3,col=ctrlcol)
      #points(xsyn,upz,type="l",lty=3,col=ctrlcol)
      if(showControls) points(xctrl,yctrl,pch=16,col=ctrlcol,cex=normcex)
      points(xvals,yvals,pch=16,col=cols,cex=cexvec)
	  points(xsel,ysel, pch=1,col=colsel,cex=cexvecsel)
	  legend("topleft",chlab, bty="n")
     #} 
}