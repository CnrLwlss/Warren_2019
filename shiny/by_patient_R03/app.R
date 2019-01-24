library(shiny)
library(data.table)
source("../plotFunctions.R")
#library(Rtsne)

library(corrgram)

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

hiliteChannel = function(dat, hilite_channel = "Z_NDUFB8"){
 hcol = colorRamp(c("red","yellow","blue"),space="Lab")
 hcol_rgb = function(x, alpha=0.3){
   vals = hcol(x)
   return(rgb(vals[1]/255,vals[2]/255,vals[3]/255,alpha))
 }
 
 dhilite = dat[dat$channel==hilite_channel,]
 dhilite$quant = ecdf(dhilite$value)(dhilite$value)
 dhcolours = sapply(dhilite$quant,hcol_rgb)
 names(dhcolours) = dhilite$cell_id
 return(dhcolours[dat$cell_id])
}

updateDat = function(dat, dtype, dsubject, dhichan, cord = c()){

    if(dtype == "2Dmito") {inpt = "Mean intensity"}else{inpt = dtype}
    dvals = dat[(dat$patrep_id==dsubject)&(dat$type==inpt),]
	
	zvals = dat[(dat$patrep_id==dsubject)&(dat$type=="z-score"),]
	
	zs = zvals$value
	names(zs) = paste(zvals$ch,zvals$cell_id)
	dvals$z = zs[paste(dvals$ch,dvals$cell_id)]	
	
    datmat = makedatmat(dvals)
	if(dhichan != " "){
	  allchans = unique(dvals$channel)
	  actual_chan = allchans[agrep(dhichan,allchans)[1]]
	  dvals$hcol = hiliteChannel(dvals, actual_chan)
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
	dvals$regression_diff = makeCond(dat,dsubject)[paste(dvals$ch,dvals$cell_id)]
	dvals$outlier_diff = factor(dvals$outlier_diff,levels=c("ABOVE","NODIFF","BELOW"))
    dvals$regression_diff = factor(dvals$regression_diff,levels=c("ABOVE","NODIFF","BELOW"))
	dvals$z_diff = factor(dvals$z_diff,levels=c("ABOVE","NODIFF","BELOW"))
	dvals
}

makeCond = function(dat,dsubject){
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
	  
#	  bdf = data.frame(patrep_id=unique(pimc$patrep_id),xvals = xall,yvals=yall)
#	  bps = brushedPoints(bdf, input$brush_main, "xvals", "yvals")

      rng = range(c(xall,yall))
      xsyn = seq(min(rng),max(rng),length.out=50)
      mod = lm(yctrl~xctrl)

      pred = predict(mod,newdata = data.frame(xctrl=xsyn), se.fit=TRUE,  interval = "prediction",na.action=na.omit)$fit
	  mid = pred[,1]
      up = pred[,3]
      low = pred[,2]
      psd = (up - low)/(2*1.96)
      upz = mid+3*psd
      lowz = mid-3*psd
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

overlaps = function(dat, dtype, dsubject, dhichan, row_condition, col_condition, cord = c()){

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
 pat = pat[,c("cell_id","channel",eval(cats))]
 patw = reshape(pat, idvar="cell_id", timevar="channel",direction="wide")
 colnames(patw) = gsub(paste(cats,".",sep=""),"",colnames(patw))
 rowcond = as.matrix(patw[,2:length(colnames(patw))]==row_condition)
 colcond = as.matrix(patw[,2:length(colnames(patw))]==col_condition)

 100*round(crossprod(rowcond,colcond)/dim(patw)[1],3)[cord,cord]

}


dat = fread("../dat.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
dat$hcol = hiliteChannel(dat)

dat$ch = substring(dat$channel,regexpr("\\_[^\\_]*$", dat$channel)+1,nchar(dat$channel))
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

dat = dat[dat$replicate==3,]

subtext =c("healthy control","nuclear-encoded mutation in CI","single, large-scale mtDNA deletion","point mutation in mito. encoded tRNA Leucine 1 (MT-TL1)","point mutation in mito. encoded tRNA (MT-TE)","point mutation in mito. encoded tRNA (MT-TG)","point mutation in mito. encoded tRNA (MT-TW)")
names(subtext) = c("Control", "CI", "Deletion", "MT-TL1", "MT-TE", "MT-TG", "MT-TW")

subdf = unique(dat[,c("patrep_id","subject_group")])
subjs = subdf$patrep_id
grps = subdf$subject_group
names(grps) = subjs
grps = grps[sort(names(grps))]

#labs = paste(subjs," (",grps[subjs],")",sep="")
labs = paste(subjs," (",subtext[grps[subjs]],")",sep="")
names(subjs) = labs
subjs = sort(subjs)


cutcords = c(2.5,3.5,4.5,6.5,7.5)#,8.5)
cordlabs = c("CI","CII","CIII","CIV","CV","OMM")#,"Cell")
cord = c("NDUFB8","GRIM19","SDHA","UqCRC2","COX4+4L2","MTCO1","OSCP","VDAC1")#,"Dystrophin","DNA1")
chlabs = c("CI","CI","CII","CIII","CIV","CIV","CV","OMM")
names(chlabs) = cord
mitochan = "VDAC1"

dat = dat[dat$ch%in%cord,]
dat$chstr = dat$ch
transform = log
dat_r = dat[dat$type=="Mean intensity",]
dat_r$type = "r (VDAC1)"
dat_theta = dat[dat$type=="Mean intensity",]
dat_theta$type = "theta (VDAC1)"

for(pid in unique(dat$patrep_id)){
 for(ch in cord){
	dt = dat[(dat$patrep_id==pid)&(dat$type=="Mean intensity"),]

	isch = as.character(dt$ch)==ch
	ismito = as.character(dt$ch)== mitochan
	prot = dt[isch,]
	mito = dt[ismito,]

	x = mito$value
	y = prot$value
	dat_r$value[(dat_r$patrep_id==pid)&(as.character(dat_r$ch)==ch)] = sqrt(x^2+y^2)
	dat_theta$value[(dat_theta$patrep_id==pid)&(as.character(dat_theta$ch)==ch)] = 360*atan(y/x)/(2*pi)
 }
}
dat=rbind(dat,dat_r,dat_theta)

types = unique(dat$type)
types = types[!types%in%c("Log mean intensity","Median intensity","Ratio median intensity (VDAC1)","Ratio log mean intensity (VDAC1)","z-score","r (VDAC1)")]
types = c("2Dmito",types)

dat$ch = factor(dat$ch, levels = cord)

jwidth = 0.25
nums = 1:length(cord) - jwidth
names(nums) = cord

dat$num = nums[dat$ch]
dat = dat[!duplicated(dat[,c("cell_id","channel","type")]),]
dat$jit = dat$num + runif(length(dat$num),0,jwidth*2) - 0.25*jwidth
dat$jitr = dat$num + runif(length(dat$num),0,jwidth) + 1.5*jwidth
dat$jitl = dat$num + runif(length(dat$num),0,jwidth) - 0.5*jwidth

nmax = 20

ui <- function(request){
fluidPage(
  titlePanel("plotIMC - interactive visualisation of imaging mass cytometry data"),
  sidebarLayout(
    sidebarPanel(
	width = 3,
	p("The default view of data from an IMC experiment is an array of interactive scatterplots.  These show protein expression levels in single skeletal fibres from one patient (coloured), compared against the same protein expression levels from all control subjects (grey).  Each panel compares the expression of a protein on the y-axis with a surrogate for mitochondrial mass on the x-axis.  Each fibre observed is represented by a point in each panel.  Individual patient fibres can be highlighted across all panels by selecting coloured points in any one panel."),
	p("Switching 'Measure of protein expression' to any option besides the default ('2Dmito') displays a stripchart, representing protein expression levels observed for a patient (coloured) compared with those observed in control subjects (grey).  Selecting patient fibres causes expression profiles for those fibres to be overlaid on top of the stripchart."),
	p("The plot below shows a matrix of Pearson's correlation coefficients between expression levels of each pair of proteins for all fibres from the selected patient.  Note that, for 2D 2Dmito plots, correlation between Ratio values is displayed."),
	p("The tables below show quantitative data and metadata summaries for all fibres as well as expression levels and biochemical deficiency categories for selected fibres (if any)."),
	p("Drag rectangles on the right panel(s) to select & highlight fibres across all channels.  To clear a selection, select an empty space on the plot."),
	p("Patient fibres are coloured according to expression in channel selected in 'Colour fibres by channel' drop-down menuL red fibres are the most highly expressed, blue fibres are the least expressed.  Note that, for 2D 2Dmito plots, colouring is by Ratio"),
	p("Data from Warren et al. (2018): Imaging mass cytometry to explore multi-dimensional respiratory chain deficiency phenotypes in single skeletal muscle fibres"),
	checkboxInput("showControls", label = "Show all control data alongside patient data?", value = TRUE, width = NULL),
    selectInput("subject", "Subject/patient", subjs),
	selectInput("type", "Measure of protein expression", types, selected="2Dmito"),#"Ratio mean intensity (VDAC1)"),
	#selectInput("nclust", "Number of k-means clusters", c(" ",as.character(1:12)), selected=""),
	#checkboxInput("showCentres", label = "Show cluster centres?", value = TRUE, width = NULL),
	#p("Only displayed if highlight channel not selected below"),
	selectInput("hichan","Colour fibres by channel",c(" ",cord), selected="NDUFB8"),
	#p("Only displayed if number of clusters not selected above.  All fibres are coloured (across all strips/channels) according to the measure observed in the selected channel"),
	fluidRow(splitLayout(cellWidths = c("33%", "34%","33%"),
	  downloadButton("download", 'Get .pdf'),
	  downloadButton("download_png", 'Get .png'),
	  bookmarkButton()
	))
    ),
	mainPanel(
   conditionalPanel(
    condition = "input.type != '2Dmito'",
    plotOutput("IMC_mainplot",height=950,
      brush = brushOpts(id = "brush_main", delay = 3000, delayType = "debounce", resetOnNew = TRUE)
   )),
   conditionalPanel(
    condition = "input.type == '2Dmito'",  
     fluidRow(
	   splitLayout(cellWidths=c("25%","25%","25%","25%"),
	   plotOutput("NDUFB8",brush = brushOpts(id = "brush_NDUFB8", delay = 3000, delayType = "debounce", resetOnNew = TRUE)),
	   plotOutput("GRIM19",brush = brushOpts(id = "brush_GRIM19", delay = 3000, delayType = "debounce", resetOnNew = TRUE)),
	   plotOutput("SDHA",brush = brushOpts(id = "brush_SDHA", delay = 3000, delayType = "debounce", resetOnNew = TRUE)),
	   plotOutput("UqCRC2",brush = brushOpts(id = "brush_UqCRC2", delay = 3000, delayType = "debounce", resetOnNew = TRUE))
	 ))
	 ),
   conditionalPanel(
    condition = "input.type == '2Dmito'",     
	 fluidRow(
	   splitLayout(cellWidths=c("25%","25%","25%","25%"),
	   plotOutput("COX4",brush = brushOpts(id = "brush_COX4", delay = 3000, delayType = "debounce", resetOnNew = TRUE)),
	   plotOutput("MTCO1",brush = brushOpts(id = "brush_MTCO1", delay = 3000, delayType = "debounce", resetOnNew = TRUE)),
	   plotOutput("OSCP",brush = brushOpts(id = "brush_OSCP", delay = 3000, delayType = "debounce", resetOnNew = TRUE))
	 ))
	 ), 
	checkboxInput("axrngCheck", label = "Allow plot y-axis range to vary between patients?", value = TRUE, width = NULL),
	htmlOutput("axrngUI"),
	h4("% fibres relative to regression 95% predictive interval"),
	tableOutput("contingency_regression"),
	h4("% fibres more likely to be from patients than controls"),
	tableOutput("contingency_outlier"),
	h4("% fibres overlapping in two categories (BELOW, NODIFF or ABOVE) from two channels"),
	fluidRow(splitLayout(cellWidths = c("25%", "25%"),
	  selectInput("overlapRows", "Overlap row:", c("BELOW","NODIFF","ABOVE"),selectize=FALSE),
      selectInput("overlapColumns", "Overlap column:", c("BELOW","NODIFF","ABOVE"),selectize=FALSE)
	  )
    ),
	tableOutput("overlap"),
	h4("Summary of fibres from selected patient:"),
	tableOutput("summdat"),
	h4("Summary of fibres from all controls:"),
	tableOutput("summdatc"),
	plotOutput("IMC_cormat",height=850),
	checkboxInput("axrngCorr", label = "Fix scatterplot axis ranges for each channel?", value = FALSE, width = NULL),
	h4("Selected fibres: summary"),
	tableOutput("summdatsel"),
	h4("Selected fibres: values"),
	tableOutput("selected_value"),
	h4("Selected fibres: regression 95% predictive interval"),
	tableOutput("selected_regression"),
	h4("Selected fibres: more likely to be from patient"),
	tableOutput("selected_outlier"),
	h4("Selected fibres: z-score cutoff"),
	tableOutput("selected_z")
	)
  )
)
}

server <- function(input, output, session) {

  onBookmark(function(state) {
    state$values$selectids <- selected$ids
  })

  # Read values from state$values when we restore bookmark
  onRestore(function(state) {
    selected$ids <- state$values$selectids
  })

plotIMC.shadecor = function (x, y, corr = NULL, col.regions, cor.method, digits = 2, cex.cor, ...) 
{
    if (is.null(corr)) {
        if (sum(complete.cases(x, y)) < 2) {
            warning("Need at least 2 complete cases for cor()")
            return()
        }
        else {
            corr <- cor(x, y, use = "pair", method = cor.method)
        }
    }
    auto <- missing(cex.cor)
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    ncol <- 140
    pal <- col.regions(ncol)
    col.ind <- as.numeric(cut(corr, breaks = seq(from = -1, to = 1, 
        length.out = ncol + 1), include.lowest = TRUE))
    abscorr <- formatC(abs(corr), digits = digits, format = "f")
    corrtxt <- formatC(corr, digits = digits, format = "f")
    if (auto) 
        cex.cor <- 0.7/strwidth(abscorr)
    usr <- par("usr")
    rect(usr[1], usr[3], usr[2], usr[4], col = pal[col.ind], 
        border = NA)
    text(0.5, 0.5, corrtxt, cex = cex.cor, col = ifelse(abs(corr)<0.5,"black","white"))
    box(col = "lightgray")
}

plotIMC.pts = function (x, y, corr = NULL, col.regions, cor.method, cex = 0.4,...) 
{
    if (!is.null(corr)) 
        return()
	cols = rgb(0,0,0,0.1)
	cols = d()$hcol
	ids = selected$ids
	hilite = d()$cell_id%in%ids
    plot.xy(xy.coords(x, y), type = "p", pch=16, col=cols,cex=cex)
	points(x[hilite],y[hilite],col="black",pch=1,lwd=0.5,cex=1.5*cex)
    box(col = "lightgray")
}

  observeEvent(input$type, {
    dat$cluster = 1
   }, priority = 1000)

  observeEvent(input$subject, {
    session$resetBrush("brush_main")
    vals$pts = dat[dat$channel=="NOCHANNEL",]
    vals$tpts = dat[dat$channel=="NOCHANNEL",]
  }, priority = 1000)
   
  vals <- reactiveValues(
	 km = NULL,
	 jitpat = "jitr",
	 jitctrl = "jitl"
  )
  
  selected <- reactiveValues(ids = c())
  
  observeEvent(eventExpr = input$brush_main, handlerExpr = {
    bps = brushedPoints(d(), input$brush_main, vals$jitpat, "value")
	if(length(bps$cell_id)==0) {
	  selected$ids=c()
	}else{
	  selected$ids = sort(unique(c(selected$ids, bps$cell_id)))
	}
  })
  
  
  observeEvent(eventExpr = input$brush_NDUFB8, handlerExpr = {
    dl = d()
	dl$value = log(dl$value)
    dw = reshape(dl[,c("cell_id","ch","value")],idvar="cell_id",timevar="ch",direction="wide")
    bps = brushedPoints(dw, input$brush_NDUFB8, paste("value",mitochan,sep="."), "value.NDUFB8")
	if(length(bps$cell_id)==0) {
	  selected$ids=c()
	}else{
	  selected$ids = sort(unique(c(selected$ids, bps$cell_id)))
	}
  })
  
    observeEvent(eventExpr = input$brush_GRIM19, handlerExpr = {
    dl = d()
	dl$value = log(dl$value)
    dw = reshape(dl[,c("cell_id","ch","value")],idvar="cell_id",timevar="ch",direction="wide")
    bps = brushedPoints(dw, input$brush_GRIM19, paste("value",mitochan,sep="."), "value.GRIM19")
	if(length(bps$cell_id)==0) {
	  selected$ids=c()
	}else{
	  selected$ids = sort(unique(c(selected$ids, bps$cell_id)))
	}
  })
  
    observeEvent(eventExpr = input$brush_SDHA, handlerExpr = {
    dl = d()
	dl$value = log(dl$value)
    dw = reshape(dl[,c("cell_id","ch","value")],idvar="cell_id",timevar="ch",direction="wide")
    bps = brushedPoints(dw, input$brush_SDHA, paste("value",mitochan,sep="."), "value.SDHA")
	if(length(bps$cell_id)==0) {
	  selected$ids=c()
	}else{
	  selected$ids = sort(unique(c(selected$ids, bps$cell_id)))
	}
  })
  
    observeEvent(eventExpr = input$brush_UqCRC2, handlerExpr = {
    dl = d()
	dl$value = log(dl$value)
    dw = reshape(dl[,c("cell_id","ch","value")],idvar="cell_id",timevar="ch",direction="wide")
    bps = brushedPoints(dw, input$brush_UqCRC2, paste("value",mitochan,sep="."), "value.UqCRC2")
	if(length(bps$cell_id)==0) {
	  selected$ids=c()
	}else{
	  selected$ids = sort(unique(c(selected$ids, bps$cell_id)))
	}
  })
  
    observeEvent(eventExpr = input$brush_COX4, handlerExpr = {
    dl = d()
	dl$value = log(dl$value)
    dw = reshape(dl[,c("cell_id","ch","value")],idvar="cell_id",timevar="ch",direction="wide")
    bps = brushedPoints(dw, input$brush_COX4, paste("value",mitochan,sep="."), "value.COX4+4L2")
	if(length(bps$cell_id)==0) {
	  selected$ids=c()
	}else{
	  selected$ids = sort(unique(c(selected$ids, bps$cell_id)))
	}
  })
  
    observeEvent(eventExpr = input$brush_MTCO1, handlerExpr = {
    dl = d()
	dl$value = log(dl$value)
    dw = reshape(dl[,c("cell_id","ch","value")],idvar="cell_id",timevar="ch",direction="wide")
    bps = brushedPoints(dw, input$brush_MTCO1, paste("value",mitochan,sep="."), "value.MTCO1")
	if(length(bps$cell_id)==0) {
	  selected$ids=c()
	}else{
	  selected$ids = sort(unique(c(selected$ids, bps$cell_id)))
	}
  })
  
    observeEvent(eventExpr = input$brush_OSCP, handlerExpr = {
    dl = d()
	dl$value = log(dl$value)
    dw = reshape(dl[,c("cell_id","ch","value")],idvar="cell_id",timevar="ch",direction="wide")
    bps = brushedPoints(dw, input$brush_OSCP, paste("value",mitochan,sep="."), "value.OSCP")
	if(length(bps$cell_id)==0) {
	  selected$ids=c()
	}else{
	  selected$ids = sort(unique(c(selected$ids, bps$cell_id)))
	}
  })
  
   observeEvent(eventExpr = input$showControls, handlerExpr = {  
	if(input$showControls) {
	  vals$jitpat = "jitr"
	  vals$jitctrl = "jitl"
	}else{
	  vals$jitpat = "jit"
	  vals$jitctrl = "jit"	  
	}
  })
  
  cond <- reactive({
    makeCond(dat,input$subject)
  })
  
  ratdat <- reactive({
    dvals = dat[(dat$patrep_id==input$subject)&(dat$type=="Ratio mean intensity (VDAC1)"),]
    if(input$hichan != " "){
	  allchans = unique(dvals$channel)
	  actual_chan = allchans[agrep(input$hichan,allchans)[1]]
	  dvals$hcol = hiliteChannel(dvals, actual_chan)
	}
	dvals
  })

  d <- reactive({
   updateDat(dat, input$type, input$subject, input$hichan, cord)
  })
  
  ctrld <- reactive({
    if(input$type == "2Dmito") {inpt = "Mean intensity"}else{inpt = input$type}
    cvals = dat[(dat$subject_group=="Control")&(dat$type==inpt),]
	cvals
  })

  output$axrngUI <- renderUI({
     if(!input$axrngCheck){
	   if(input$type == "2Dmito") {inpt = "Mean intensity"}else{inpt = input$type}
	   dsub = dat[dat$type==inpt,]
	   if(input$type == "2Dmito"){
	     values = log(dsub$value[dsub$value>0])
		 values = values[is.finite(values)]
	     minval = min(values,na.rm=TRUE)
		 maxval =max(values,na.rm=TRUE)
	   }else{
	   	 minval = min(dsub$value,na.rm=TRUE)
		 maxval = max(dsub$value,na.rm=TRUE)
	   }
	   stepval = (maxval - minval)/1000
       sliderInput("axrng", "Fix range for axes of all plots:", min = minval, max = maxval, step=stepval, width="850px", value = c(minval, maxval), dragRange = TRUE)
	  }
  })
 
  output$plot_brushedpoints = renderTable(within(d()[d()$cell_id%in%selected$ids,],rm(jitl,jitr,jit,num,chstr,hcol,colour)))
  
  output$summdat = renderTable({summtab(d(),cord)},rownames=TRUE)
  output$summdatc = renderTable({summtab(ctrld(),cord)},rownames=TRUE)
  output$summdatsel = renderTable({summtab(d()[d()$cell_id%in%selected$ids,],cord)},rownames=TRUE)
  
  output$selected_value = renderTable(getwide(d()[d()$cell_id%in%selected$ids,],"value"))
  output$selected_outlier = renderTable(getwide(d()[d()$cell_id%in%selected$ids,],"outlier_diff"))
  output$selected_regression = renderTable(getwide(d()[d()$cell_id%in%selected$ids,],"regression_diff"))
  output$selected_z = renderTable(getwide(d()[d()$cell_id%in%selected$ids,],"z_diff"))
  
  output$contingency_outlier = renderTable(as.data.frame.matrix(100*with(d(),table(outlier_diff,ch))/length(unique(d()$cell_id))),include.rownames=TRUE)
  output$contingency_regression = renderTable(as.data.frame.matrix(100*with(d(),table(regression_diff,ch))/length(unique(d()$cell_id))),include.rownames=TRUE)
  output$contingency_z = renderTable(as.data.frame.matrix(100*with(d(),table(z_diff,ch))/length(unique(d()$cell_id))),include.rownames=TRUE)
  
  output$overlap = renderTable(overlaps(dat, input$type, input$subject, input$hichan, input$overlapRows, input$overlapColumns, cord),include.rownames=TRUE)
  

	mystripchart = function(){
	 schart(d(),ctrld(),selected$ids,subtext=subtext, nums = nums, cord = cord, cutcords = cutcords, cordlabs = cordlabs,subjectLabel=input$subject,hichan=input$hichan,axrngCheck = input$axrngCheck, axrng=input$axrng,jitpat=vals$jitpat,jitctrl=vals$jitctrl, showControls = input$showControls)
	}
	
	myarrayplot = function(){
	 op = par(mfrow=c(4,2),mar=c(3.2, 3.2, 1, 0.5) + 0.1,mgp=c(1.75, 0.7, 0))
	  alld = d()
	  allc = ctrld()
	  allrat = ratdat()
	  for(i in seq_along(cord[cord!=mitochan])){
	   ch = cord[cord!=mitochan][i]
	   chlab = paste(LETTERS[i]," (",chlabs[i],")",sep="")
	   arrayplot(alld,allc,allrat,cord=cord,ch=ch,ids=selected$ids,reg_diff = cond(),mitochan="VDAC1",hichan=input$hichan,showControls = input$showControls,axrngCheck = input$axrngCheck, axrng=input$axrng,chlab=chlab)
	  }
	 par(op) 
	}
	
	mycormat = function(){
	 if(input$type == "2Dmito") {inpt = "Ratio mean intensity (VDAC1)"}else{inpt = input$type}
	 dvs = dat[(dat$type==inpt),]
     mins = aggregate(dvs$value,by=list(dvs$ch),FUN=min)
     maxs = aggregate(dvs$value,by=list(dvs$ch),FUN=max)
     diffs = maxs
     diffs$x = diffs$x - mins$x

     dscl = d()
	 if((inpt=="z-score")|(grepl(mitochan,inpt))) dscl = dscl[dscl$ch!=mitochan,]
     for(ch in mins$Group.1){
       minval = mins$x[mins$Group.1==ch]
	   diffsval = diffs$x[diffs$Group.1==ch]
       dscl$value[dscl$ch==ch] = (dscl$value[dscl$ch==ch] - minval)/diffsval
     }
	
     widef = reshape(subset(dscl,select=c("value","cell_id","ch")),idvar="cell_id",timevar="ch",v.names="value",direction="wide")
     rownames(widef)=widef$cell_id
     widef$cell_id=NULL
     colnames(widef)=gsub("value.","",colnames(widef))
     colnames(widef)=gsub("MED_","",colnames(widef))
     colnames(widef)=gsub("R_LOG_","",colnames(widef))
     colnames(widef)=gsub("LOG_","",colnames(widef))
     colnames(widef)=gsub("R_","",colnames(widef))
     colnames(widef)=gsub("Z_","",colnames(widef))
     widefnames = colnames(widef)
     cordnew = cord[cord%in%widefnames]
     widef = subset(widef,select = cordnew)
	 if(!input$axrngCorr){
	   axrng = range(dscl$value[is.finite(dscl$value)],na.rm=TRUE)
	 }else{
	   axrng = c(0,1)
	 }
     corrgram(widef,order=FALSE,lower.panel=plotIMC.shadecor,upper.panel=plotIMC.pts,abs=TRUE,col.regions = colorRampPalette(c("blue","white","red")),xlim=axrng,ylim=axrng,main=inpt)	
	}
  
  output$IMC_mainplot = 
    renderPlot({
	  mystripchart()
    },width=950,height=850,pointsize=26)
  
  #output$IMC_arrayplot <- renderPlot({
  #  myarrayplot()
  #},width=850,height=850,pointsize=20)
  
    output$NDUFB8 <- renderPlot({
	op = par(mar=c(4, 2.75, 1, 0.5) + 0.1,mgp=c(1.75, 0.7, 0))
	    arrayplot(d(),ctrld(),ratdat(),cord=cord,ch="NDUFB8",ids=selected$ids,reg_diff = cond(),mitochan="VDAC1",hichan=input$hichan,showControls=input$showControls,axrngCheck = input$axrngCheck, axrng=input$axrng, chlab=paste("A (",chlabs["NDUFB8"],")",sep=""),logify = TRUE)
	par(op)
	})
	output$GRIM19 <- renderPlot({
	op = par(mar=c(4, 2.75, 1, 0.5) + 0.1,mgp=c(1.75, 0.7, 0))
	    arrayplot(d(),ctrld(),ratdat(),cord=cord,ch="GRIM19",ids=selected$ids,reg_diff = cond(),mitochan="VDAC1",hichan=input$hichan,showControls=input$showControls,axrngCheck = input$axrngCheck, axrng=input$axrng, chlab=paste("B (",chlabs["GRIM19"],")",sep=""),logify = TRUE)
	par(op)
	})
	output$SDHA <- renderPlot({
	op = par(mar=c(4, 2.75, 1, 0.5) + 0.1,mgp=c(1.75, 0.7, 0))
	    arrayplot(d(),ctrld(),ratdat(),cord=cord,ch="SDHA",ids=selected$ids,reg_diff = cond(),mitochan="VDAC1",hichan=input$hichan,showControls=input$showControls,axrngCheck = input$axrngCheck, axrng=input$axrng, chlab=paste("C (",chlabs["SDHA"],")",sep=""),logify = TRUE)
	par(op)
	})
	output$UqCRC2 <- renderPlot({
	op = par(mar=c(4, 2.75, 1, 0.5) + 0.1,mgp=c(1.75, 0.7, 0))
	    arrayplot(d(),ctrld(),ratdat(),cord=cord,ch="UqCRC2",ids=selected$ids,reg_diff = cond(),mitochan="VDAC1",hichan=input$hichan,showControls=input$showControls,axrngCheck = input$axrngCheck, axrng=input$axrng, chlab=paste("D (",chlabs["UqCRC2"],")",sep=""),logify = TRUE)
	par(op)
	})
	output$COX4 <- renderPlot({
	op = par(mar=c(4, 2.75, 1, 0.5) + 0.1,mgp=c(1.75, 0.7, 0))
	    arrayplot(d(),ctrld(),ratdat(),cord=cord,ch="COX4+4L2",ids=selected$ids,reg_diff = cond(),mitochan="VDAC1",hichan=input$hichan,showControls=input$showControls,axrngCheck = input$axrngCheck, axrng=input$axrng, chlab=paste("E (",chlabs["COX4+4L2"],")",sep=""),logify = TRUE)
	par(op)
	})
	output$MTCO1 <- renderPlot({
	op = par(mar=c(4, 2.75, 1, 0.5) + 0.1,mgp=c(1.75, 0.7, 0))
	    arrayplot(d(),ctrld(),ratdat(),cord=cord,ch="MTCO1",ids=selected$ids,reg_diff = cond(),mitochan="VDAC1",hichan=input$hichan,showControls=input$showControls,axrngCheck = input$axrngCheck, axrng=input$axrng, chlab=paste("F (",chlabs["MTCO1"],")",sep=""),logify = TRUE)
	par(op)
	})
	output$OSCP <- renderPlot({
	op = par(mar=c(4, 2.75, 1, 0.5) + 0.1,mgp=c(1.75, 0.7, 0))
	    arrayplot(d(),ctrld(),ratdat(),cord=cord,ch="OSCP",ids=selected$ids,reg_diff = cond(),mitochan="VDAC1",hichan=input$hichan,showControls=input$showControls,axrngCheck = input$axrngCheck, axrng=input$axrng, chlab=paste("G (",chlabs["OSCP"],")",sep=""),logify = TRUE)
	par(op)
	})
  
  output$IMC_stripchart <- renderPlot({
    mystripchart()
  },width=850,height=850,pointsize=26)
  
  output$IMC_cormat <-renderPlot({
    mycormat()
  },width=850,height=850,pointsize=26)
  
   output$download <- downloadHandler(
    filename =  paste("stripchart","pdf",sep="."),
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
	  if(input$type=="2Dmito"){
	  pdf(file, pointsize=16,width = 8.27, height = 16.54)
        myarrayplot()
	  }else{
	  pdf(file, pointsize=16,width = 8.27, height = 8.27)
	    mystripchart()
	  }
      dev.off()
    }
  )	

    output$download_png <- downloadHandler(
	
    filename =  paste("stripchart","png",sep="."),
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
	  if(input$type=="2Dmito"){
       png(file, width=1750, height=3500, pointsize=51)
        myarrayplot()
	  }else{
       png(file, width=1750, height=1750, pointsize=51)
	    mystripchart()
	  }
      dev.off()  # turn the device off
    }
  )	


}

shinyApp(ui,server,enableBookmarking = "url")
