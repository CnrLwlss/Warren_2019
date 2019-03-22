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

# Make summary row for d() object
summdfunc = function(d,cord,transdat=identity,FUN,...) {
  summ =aggregate(transdat(d$value),by=list(factor(d$ch, levels = cord)),FUN=FUN,...)$x
  names(summ) = cord
  return(summ)
}

# Make summary table for d()
summtab = function(d,cord,transdat=identity){
    if(dim(d)[1]>0){
    res=list()
    res$min = summdfunc(d,cord,transdat=transdat,min)
    res[["2.5%ile"]] = summdfunc(d,cord,transdat=transdat,quantile,0.025)
    res[["5%ile"]] = summdfunc(d,cord,transdat=transdat,quantile,0.05)
    res[["25%ile"]] = summdfunc(d,cord,transdat=transdat,quantile,0.25)
    res$median = summdfunc(d,cord,transdat=transdat,median)
    res$mean = summdfunc(d,cord,transdat=transdat,mean)
    res[["75%ile"]] = summdfunc(d,cord,transdat=transdat,quantile,0.75)
    res[["95%ile"]] = summdfunc(d,cord,transdat=transdat,quantile,0.95)
    res[["97.5%ile"]] = summdfunc(d,cord,transdat=transdat,quantile,0.975)
    res$max = summdfunc(d,cord,transdat=transdat,max)
    res$SD = summdfunc(d,cord,transdat=transdat,sd)
    do.call("rbind",res)
    }else{
    data.frame()
    }
}

makedatmat = function(d){
  datw = reshape(d, idvar="cell_id",timevar="channel",direction="wide",drop=c("id","patient_id","patrep_id","colour","patient_type","subject_group","jit","cluster","type","num","ch","hcol"))
  colnames(datw) = gsub("value.","",colnames(datw))
  datmat = as.matrix(datw[,2:length(colnames(datw))])
  rownames(datmat) = datw$cell_id
  return(datmat)
}

tSNE_data = function(d,nclust=1,perp=25,theta=0.5){
  set.seed(42)

  datmat = makedatmat(d)
  tmod = Rtsne(datmat, check_duplicates = FALSE,pca = TRUE, perplexity = perp, theta = theta, dims = 2)
  dts = as.data.frame(tmod$Y, stringsAsFactors = FALSE)
  colnames(dts) = c("tSNE1","tSNE2")
  dts$cell_id = rownames(datmat)

  if(nclust>1){
    km = kmeans(scale(datmat),nclust)
    dts$cluster = km$cluster
  }else{
    dts$cluster = 1
  } 
  return(dts)
}

tSNE_plot = function(dts, cols){
  main = paste("N =",length(unique(dts$cell_id)))
  op = par(mar=c(6, 4, 2, 0) + 0.1)
  plot(dts$tSNE1,dts$tSNE2,xlab="t-SNE 1",ylab="t-SNE 2",main=main,type="n")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(0.65,0.65,0.65))
  points(dts$tSNE1,dts$tSNE2,pch=16,cex=1.0,col=cols[dts$cell_id])
  par(op)
  return(dts)
}

hiliteChannel = function(dat, hilite_ch = "NDUFB8", hilite_type = "theta (VDAC1)"){
 hcol = colorRamp(c("red","yellow","blue"),space="Lab")
 hcol_rgb = function(x, alpha=0.3){
   vals = hcol(x)
   return(rgb(vals[1]/255,vals[2]/255,vals[3]/255,alpha))
 }

 dhilite = dat[(as.character(dat$ch)==hilite_ch)&(dat$type==ifelse(hilite_type=="2Dmito","theta (VDAC1)",hilite_type)),]
 dhilite$quant = ecdf(dhilite$value)(dhilite$value)
 
 dhcolours = sapply(dhilite$quant,hcol_rgb)
 names(dhcolours) = dhilite$cell_id
 return(dhcolours[dat$cell_id])
}

updateDat = function(dat, dtype, dsubject, dhichan, cord = c(), clevel = 0.95){

    if(dtype == "2Dmito") {inpt = "Mean intensity"}else{inpt = dtype}
    dvals = dat[(dat$patrep_id==dsubject)&(dat$type==inpt),]
	
	zvals = dat[(dat$patrep_id==dsubject)&(dat$type=="z-score"),]
	
	zs = zvals$value
	names(zs) = paste(zvals$ch,zvals$cell_id)
	dvals$z = zs[paste(dvals$ch,dvals$cell_id)]	
	
    datmat = makedatmat(dvals)
	if(dhichan != " "){
	  hitype = ifelse(dtype=="2Dmito","theta (VDAC1)",dtype)
	  dvals$hcol = hiliteChannel(dat[(dat$patrep_id==dsubject)&(dat$type==hitype),], dhichan,hitype)
	  if((dtype == "z-score")||(grepl("Ratio",dtype))) dvals$hcol[dvals$ch==mitochan] = rgb(0/255,154/255,73/255,0.25)
	  if(dtype == "z-score") dvals = dvals[!dvals$ch%in%c("Dystrophin","DNA1"),]
	}
	csig = dat[(dat$subject_group=="Control")&(dat$type==inpt),]
    dvals$outlier_diff = "NODIFF"
	for(j in seq_along(cord)){
        pat = dvals$value[as.character(dvals$ch)==as.character(cord[j])]
        ctr = csig$value[as.character(csig$ch)==as.character(cord[j])]
	    diffs = findDiffs(pat,ctr)
		dvals$outlier_diff[as.character(dvals$ch)==as.character(cord[j])] = diffs
	}
	dvals$z_diff = "NODIFF"
	dvals$z_diff[dvals$z>3] = "ABOVE"
	dvals$z_diff[dvals$z<(-3)] = "BELOW"
	dvals$regression_diff = makeCond(dat, dsubject, clevel)[paste(dvals$ch,dvals$cell_id)]
	dvals$outlier_diff = factor(dvals$outlier_diff,levels=c("ABOVE","NODIFF","BELOW"))
    dvals$regression_diff = factor(dvals$regression_diff,levels=c("ABOVE","NODIFF","BELOW"))
	dvals$z_diff = factor(dvals$z_diff,levels=c("ABOVE","NODIFF","BELOW"))
	dvals
}

makeCond = function(dat,dsubject, clevel = 0.95){
	pimc = dat[(dat$patrep_id==dsubject)&(dat$type=="Mean intensity"),]
	cimc = dat[(dat$subject_group=="Control")&(dat$type=="Mean intensity"),]
	pimc$regression_diff = "NODIFF"
	
	transform = log
	for(ch in cord[cord!=mitochan]){
      xctrl = transform(cimc$value[(as.character(cimc$ch)==mitochan)&(cimc$subject_group=="Control")])
      yctrl = transform(cimc$value[(as.character(cimc$ch)==as.character(ch))&(cimc$subject_group=="Control")])
	  xvals = transform(pimc$value[as.character(pimc$ch)==mitochan])
      yvals = transform(pimc$value[as.character(pimc$ch)==ch])
	  xall = c(xctrl,xvals)
	  yall = c(yctrl,yvals)
	  N = length(yvals)
	  regression_diff = rep("NODIFF",N)

      rng = range(c(xall,yall))
      xsyn = seq(min(rng),max(rng),length.out=50)
      mod = lm(yctrl~xctrl)

      pred = predict(mod,newdata = data.frame(xctrl=xsyn), se.fit=TRUE,  interval = "prediction",na.action=na.omit, level = clevel)$fit
	  mid = pred[,1]
      up = pred[,3]
      low = pred[,2]
	  upy = approxfun(xsyn,up)
      lowy = approxfun(xsyn,low)
      
      below = yvals<lowy(xvals)
      above = yvals>upy(xvals)
	  regression_diff[below]="BELOW"
	  regression_diff[above]="ABOVE"
	  pimc$regression_diff[as.character(pimc$ch)==ch] = regression_diff
	}
	rd = pimc$regression_diff
	names(rd) = paste(pimc$ch,pimc$cell_id)
	rd
}

overlaps = function(dat, dtype, dsubject, dhichan, row_condition, col_condition, cord = c(),counts = FALSE){

 if(dtype == "theta (VDAC1)"){
   cats = "outlier_diff"
 }else if(dtype == "r (VDAC1)"){
   cats = "outlier_diff"
 }else if(dtype == "2Dmito"){
   cats = "regression_diff"
 }else if(dtype == "z-score"){
   cats = "z_diff"
 }else{
   cats = ""
 }
 
 pat = data.frame(updateDat(dat, dtype, dsubject, dhichan, cord))
 pat = pat[,c("cell_id","ch",eval(cats))]
 patw = reshape(pat, idvar="cell_id", timevar="ch",direction="wide")
 colnames(patw) = gsub(paste(cats,".",sep=""),"",colnames(patw))
 rowcond = as.matrix(patw[,2:length(colnames(patw))]==row_condition)
 colcond = as.matrix(patw[,2:length(colnames(patw))]==col_condition)
 
 rawcounts = crossprod(rowcond,colcond)[cord,cord]
 percentages = 100*round(rawcounts/dim(patw)[1],3)
 if(counts){
   return(as.data.frame.matrix(rawcounts))
 }else{
   return(as.data.frame.matrix(percentages))
 }
}

contig_outlier = function(dat,counts=FALSE){
  dfm = as.data.frame.matrix(with(dat,table(outlier_diff,ch)))
  if(counts){return(dfm)}else{return(100*dfm/length(unique(dat$cell_id)))}
}

contig_regression = function(dat,counts=FALSE){
  dfm = as.data.frame.matrix(with(dat,table(regression_diff,ch)))
  if(counts){return(dfm)}else{return(100*dfm/length(unique(dat$cell_id)))}
}

contig_z = function(dat,counts=FALSE){
  dfm = as.data.frame.matrix(with(dat,table(z_diff,ch)))
  if(counts){return(dfm)}else{return(100*dfm/length(unique(dat$cell_id)))}
}


