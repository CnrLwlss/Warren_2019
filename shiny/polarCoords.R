head(dat)

# Find values which are more likely to be from patient than control
findDiffs = function(dv,cv,low=0.05,high=0.95){
 res = rep("NODIFF",length(dv))
 allv = c(dv,cv)
 allv = allv[!is.na(allv)]
 if(length(unique(allv))>1){
  dvd = approxfun(density(dv,from=min(allv),to=max(allv)))
  cvd = approxfun(density(cv,from=min(allv),to=max(allv)))
  lower = quantile(cv,low)
  upper = quantile(cv,high)
 
  res[(dvd(dv)>cvd(dv))&(dv>upper)]="ABOVE"
  res[(dvd(dv)>cvd(dv))&(dv<lower)]="BELOW"
 }
 return(res)
}

cord = c("NDUFB8","GRIM19","SDHA","UqCRC2","COX4+4L2","MTCO1","OSCP","VDAC1")#,"Dystrophin","DNA1")
mitochan = "VDAC1"

pids = sort(unique(dat$patient_id))

lowval=0.025
highval=0.975
transform = log

pdf("Polar.pdf")

for(pid in pids){
 for(ch in cord[cord!=mitochan]){
	dt = dat[(dat$patient_id==pid)&(dat$type=="Mean intensity"),]
	cdt = dat[(dat$patient_type=="control")&(dat$type=="Mean intensity"),]
	rng = range(transform(dat$value[(dat$type=="Mean intensity")&(as.character(dat$ch)==ch)]),na.rm=TRUE)
      mrng = range(transform(dat$value[(dat$type=="Mean intensity")&(as.character(dat$ch)==mitochan)]),na.rm=TRUE)

	isch = as.character(dt$ch)==ch
	ismito = as.character(dt$ch)== mitochan
	prot = dt[isch,]
	mito = dt[ismito,]
	
	cisch = as.character(cdt$ch)==ch
	cismito = as.character(cdt$ch)== mitochan
	cprot = cdt[cisch,]
	cmito = cdt[cismito,]

	x = transform(mito$value)
	y = transform(prot$value)
	r = sqrt(x^2+y^2)
	theta = atan(y/x)


	cx = transform(cmito$value)
	cy = transform(cprot$value)
	cr = sqrt(cx^2+cy^2)
	ctheta = atan(cy/cx)

	
	coords = data.frame(x=x,y=y,theta=theta,r=r)
	ccoords = data.frame(x=cx,y=cy,theta=ctheta,r=cr)

	dens = density(coords$theta)
	cdens = density(ccoords$theta)

	coords$diffs = findDiffs(coords$theta,ccoords$theta,low=lowval,high=highval)
	N = length(coords$diffs)
	above = sum(coords$diffs=="ABOVE")
	below = sum(coords$diffs=="BELOW")
	coords$col = rgb(1,0,0,0.2)
	coords$col[coords$diff%in%c("ABOVE","BELOW")] = rgb(0,0,1,0.2)

	mlab = paste(pid,ch,"\nN:",N,"Above:",above,"Below:",below)

	op=par(mfrow=c(2,2))
	 plot(ccoords$x,ccoords$y,xlab=mitochan,ylab=ch,main=paste("Cartesian",mlab),col=rgb(0,0,0,0.2),pch=16,xlim=mrng,ylim=rng,cex=0.5)
       abline(a=0,b=tan(quantile(ccoords$theta,lowval)),lty=2,col=rgb(0,0,0,0.2),lwd=2)
       abline(a=0,b=tan(quantile(ccoords$theta,highval)),lty=2,col=rgb(0,0,0,0.2),lwd=2)
	 abline(a=0,b=tan(quantile(ccoords$theta,0.5)),lty=1,col=rgb(0,0,0,0.2),lwd=2)
	 points(coords$x,coords$y,col=coords$col,pch=16,cex=0.75)

	 trng = range(c(coords$theta,ccoords$theta))
 	 rrng = range(c(coords$r,ccoords$r))  
	 plot(ccoords$theta,ccoords$r,xlab="theta (rad)",ylab="r",main=paste("Polar",mlab),col=rgb(0,0,0,0.2),pch=16,xlim=trng,ylim=rrng,cex=0.5)
	 abline(v=quantile(ccoords$theta,c(lowval,highval)),lty=2,col=rgb(0,0,0,0.2),lwd=2)
	 points(coords$theta,coords$r,col=coords$col,pch=16,cex=0.75)	 

	 plot(dens,xlab = "theta (rad)",main="Density",xlim=range(c(dens$x,cdens$x)),ylim=range(c(dens$y,cdens$y)),col=rgb(1,0,0,0.6),lwd=2)
	 points(cdens,type="l",col=rgb(0,0,0,0.6),lwd=2)
	 pts = coords$theta[coords$diff%in%c("ABOVE","BELOW")]
	 points(pts,rep(0,length(pts)),col= rgb(0,0,1,0.2),pch=3,lwd=2)
	 abline(v=quantile(ccoords$theta,c(lowval,highval)),lty=2,col=rgb(0,0,0,0.6),lwd=2)
	 legend("topleft",c("Controls","Patient"),col=c(rgb(0,0,0,0.2),rgb(1,0,0,0.6)),bg="white",cex=0.75,lwd=1)

	 densr = density(coords$r)
	 cdensr = density(ccoords$r)

	 plot(densr,xlab = "r",main="Density",xlim=range(c(densr$x,cdensr$x)),ylim=range(c(densr$y,cdensr$y)),col=rgb(1,0,0,0.6),lwd=2)
	 points(cdensr,type="l",col=rgb(0,0,0,0.6),lwd=2)
	 abline(v=quantile(ccoords$r,c(lowval,highval)),lty=2,col=rgb(0,0,0,0.6),lwd=2)
	 legend("topleft",c("Controls","Patient"),col=c(rgb(0,0,0,0.2),rgb(1,0,0,0.6)),bg="white",cex=0.75,lwd=1)

	par(op)
 }
}
dev.off()
 
