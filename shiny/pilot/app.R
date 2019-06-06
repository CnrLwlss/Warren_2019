library(shiny)
library(data.table)
library(corrgram)

source("../plotFunctions.R", local = TRUE)
source("../dataFunctions.R", local = TRUE)

subtext =c("healthy control","POLG")
names(subtext) = c("Control", "Patient")

cutcords = c(2.5,3.5,4.5,6.5,7.5)#,8.5)
cordlabs = c("CI","CII","CIII","CIV","CV","OMM")#,"Cell")
cord = c("NDUFA13","NDUFB8","SDHA","UqCRC2","COX4+4L2","MTCO1","OSCP","VDAC1")#,"Dystrophin","DNA1")
chlabs = c("CI","CI","CII","CIII","CIV","CIV","CV","OMM")
names(chlabs) = cord
mitochan = "VDAC1"

source("../parseData.R", local = TRUE)

fulldat = "../pilot_data.txt"
zipdat = "../pilot_datasmall.txt.gz"

cord = c("NDUFA13","NDUFB8","SDHA", "UqCRC2","COX4+4L2","MTCO1","OSCP","VDAC1")#,"TOM20","Dystrophin","DNA1","DNA2")
chlabs=c("CI","CI","CII", "CIII","CIV","CIV","CV","OMM")#,"OMM","OMM","OMM","OMM")
names(chlabs) = cord
mitochan = "VDAC1"

if(!file.exists(zipdat)){
 dat = getData(fulldat,cord,mitochan)
 write.table(dat,gzfile(zipdat),quote=FALSE,sep="\t",row.names=FALSE)
}else{
 dat = fread(zipdat,sep="\t",stringsAsFactors=FALSE,header=TRUE)
}

if(grepl("R03",basename(getwd()))){
  repnum = 3
}else if(grepl("R02",basename(getwd()))){
  repnum = 2
}else{
  repnum = 1
}

dat = dat[dat$replicate==repnum,]

subdf = unique(dat[,c("patrep_id","subject_group")])
subjs = subdf$patrep_id
grps = subdf$subject_group
names(grps) = subjs
grps = grps[sort(names(grps))]

#labs = paste(subjs," (",grps[subjs],")",sep="")
labs = paste(subjs," (",subtext[grps[subjs]],")",sep="")
names(subjs) = labs
subjs = sort(subjs)

dat$hcol = hiliteChannel(dat,hilite_ch = "NDUFB8",hilite_type = "theta (VDAC1)")

types = unique(dat$type)
types = types[!types%in%c("Log mean intensity","Median intensity","Ratio median intensity (VDAC1)","Ratio mean intensity (VDAC1)","Ratio log mean intensity (VDAC1)","z-score","r (VDAC1)")]
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
	p("This webpage is a tool for interactive analysis of IMC data gathered from skeletal muscle fibre sections sampled from patients with mitochondrial disease.  Pilot data for Warren et al. (2019): An Imaging Mass Cytometric approach to decoding mitochondrial heterogeneity at a single muscle fibre resolution"),
	checkboxInput("showControls", label = "Show all control data alongside patient data?", value = TRUE, width = NULL),
    selectInput("subject", "Subject/patient", subjs),
	selectInput("type", "Measure of protein expression", types, selected="2Dmito"),#"Ratio mean intensity (VDAC1)"),
	selectInput("hichan","Colour fibres by channel",c(" ",cord), selected="NDUFB8"),
	fluidRow(splitLayout(cellWidths = c("33%", "34%","33%"),
	  downloadButton("download", 'Get .pdf'),
	  downloadButton("download_png", 'Get .png'),
	  bookmarkButton()
	)),
	p(""),
	p("The default view of data is an array of interactive scatterplots: 2Dmito, comparing the expression of a protein on the y-axis with a surrogate for mitochondrial mass on the x-axis.  Coloured points represent fibres from a patient, and grey points represent fibres from all control subjects.  Each fibre observed is represented by a point in each scatterplot.  Solid grey line is linear regression through control data.  Dashed lines represent boundaries of 95% predcitive interval for control fibres.  Points are coloured according to expression of the proteins selected in the 'Colour fibres by channel' drop-down menu: red fibres express the selected protein least, blue fibres express that protein the most.  Amount of expression is measured by theta (see below).  To emphasise the multiple measurements made for each fibre, selecting fibres on any one plot (including stripcharts) causes circles to be drawn around the position of those fibres in all other plots."),
	p("Switching 'Measure of protein expression' to any option besides the default (2Dmito) displays a stripchart, representing the distributions of protein expression levels observed for the selected patient (coloured) compared with those observed in control subjects (grey)."),
	p("To select & highlight the expression of all proteins for selected fibres, drag rectangles on the plots using the mouse.  To clear a selection, select any empty space on the plot."),
	p("Two tables below the main panel summarise the proportion of fibres belonging to each of three categories: sigificantly ABOVE, BELOW or not different from (NODIFF) control fibres, for each protein.  A third table summarises all two-way combinations of overlap between channels for any pair of the three categories listed above (select category combinations using drop-down menus).  Two more tables summarise expression levels for each protein: one for the selected patient and another for all controls."),
	p("The panel below shows a matrix of Pearson's correlation coefficients between expression levels (theta) of each pair of proteins for all fibres from the selected patient."),
	p("Below that, there are four tables summarising the properties of any selected fibres.")
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
	   plotOutput("NDUFA13",brush = brushOpts(id = "brush_NDUFA13", delay = 3000, delayType = "debounce", resetOnNew = TRUE)),
	   plotOutput("NDUFB8",brush = brushOpts(id = "brush_NDUFB8", delay = 3000, delayType = "debounce", resetOnNew = TRUE)),
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
	checkboxInput("counts", label = "Show raw count data in table (alternative is percentages)?", value = FALSE, width = NULL),
	htmlOutput("axrngUI"),
	h4("Proportion (%) of fibres lying outside predictive interval for control fibres in 2Dmito plot"),
	tableOutput("contingency_regression"),
	h4("Proportion (%) of fibres more likely to be from patients than controls (measured using theta)"),
	tableOutput("contingency_outlier"),
	h4("Proportion (%) of fibres overlapping between selected categories (BELOW, NODIFF or ABOVE) from all pairwise combinations of proteins (rows & columns)"),
	fluidRow(splitLayout(cellWidths = c("25%", "25%"),
	  selectInput("overlapRows", "Overlap row:", c("BELOW","NODIFF","ABOVE"),selectize=FALSE),
      selectInput("overlapColumns", "Overlap column:", c("BELOW","NODIFF","ABOVE"),selectize=FALSE)
	  )
    ),
	tableOutput("overlap"),
	h4("Summary of distribution of expression levels for fibres from selected patient:"),
	tableOutput("summdat"),
	h4("Summary of distribution of expression levels for fibres from all controls:"),
	tableOutput("summdatc"),
	plotOutput("IMC_cormat",height=850),
	checkboxInput("axrngCorr", label = "Fix scatterplot axis ranges for each channel?", value = FALSE, width = NULL),
	h4("Selected fibres: summary"),
	tableOutput("summdatsel"),
	h4("Selected fibres: raw values"),
	tableOutput("selected_value"),
	h4("Selected fibres: fibre categories, predictive interval for contols"),
	tableOutput("selected_regression"),
	h4("Selected fibres: fibre categories, (theta) more likely to be from patients than controls"),
	tableOutput("selected_outlier")
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
	#points(x[hilite],y[hilite],col="black",pch=1,lwd=0.5,cex=1.5*cex)
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

  observeEvent(eventExpr = input$brush_NDUFA13, handlerExpr = {
    dl = d()
	dl$value = log(dl$value)
    dw = reshape(dl[,c("cell_id","ch","value")],idvar="cell_id",timevar="ch",direction="wide")
    bps = brushedPoints(dw, input$brush_NDUFA13, paste("value",mitochan,sep="."), "value.NDUFA13")
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
    #dvals = dat[(dat$patrep_id==input$subject)&(dat$type=="Ratio mean intensity (VDAC1)"),]
    dvals = dat[(dat$patrep_id==input$subject)&(dat$type=="theta (VDAC1)"),]
    if(input$hichan != " "){
	  dvals$hcol = hiliteChannel(dvals, input$hichan, "theta (VDAC1)")
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
  
  makeOverlaps <- reactive({overlaps(dat, ifelse(input$type=="Mean intensity","2Dmito",input$type), input$subject, input$hichan, input$overlapRows, input$overlapColumns, cord, input$counts)})
  makeContOut <- reactive({contig_outlier(d(),input$counts)})
  makeContReg <- reactive({contig_regression(d(),input$counts)})
  makeContZ <- reactive({contig_z(d(),input$counts)})
 
  output$plot_brushedpoints = renderTable(within(d()[d()$cell_id%in%selected$ids,],rm(jitl,jitr,jit,num,chstr,hcol,colour)))
  
  output$summdat = renderTable({summtab(d(),cord)},rownames=TRUE)
  output$summdatc = renderTable({summtab(ctrld(),cord)},rownames=TRUE)
  output$summdatsel = renderTable({summtab(d()[d()$cell_id%in%selected$ids,],cord)},rownames=TRUE)
  
  output$selected_value = renderTable(getwide(d()[d()$cell_id%in%selected$ids,],"value"))
  output$selected_outlier = renderTable(getwide(d()[d()$cell_id%in%selected$ids,],"outlier_diff"))
  output$selected_regression = renderTable(getwide(d()[d()$cell_id%in%selected$ids,],"regression_diff"))
  output$selected_z = renderTable(getwide(d()[d()$cell_id%in%selected$ids,],"z_diff"))
  
  #output$contingency_outlier = renderTable(as.data.frame.matrix(100*with(d(),table(outlier_diff,ch))/length(unique(d()$cell_id))),include.rownames=TRUE)
  #output$contingency_regression = renderTable(as.data.frame.matrix(100*with(d(),table(regression_diff,ch))/length(unique(d()$cell_id))),include.rownames=TRUE)
  #output$contingency_z = renderTable(as.data.frame.matrix(100*with(d(),table(z_diff,ch))/length(unique(d()$cell_id))),include.rownames=TRUE)
  
  output$contingency_outlier = renderTable(makeContOut(),include.rownames=TRUE)
  output$contingency_regression = renderTable(makeContReg(),include.rownames=TRUE)
  output$contingency_z = renderTable(makeContZ(),include.rownames=TRUE)
  
  output$overlap = renderTable(makeOverlaps(),include.rownames=TRUE)
  

	mystripchart = function(){
	 schart(d(),ctrld(),selected$ids,subtext=subtext, nums = nums, cord = cord, cutcords = cutcords, cordlabs = cordlabs,subjectLabel=input$subject,hichan=input$hichan,axrngCheck = input$axrngCheck, axrng=input$axrng,jitpat=vals$jitpat,jitctrl=vals$jitctrl, showControls = input$showControls)
	}
	
	myarrayplot = function(){
	 op = par(mfrow=c(4,3),mar=c(3.2, 3.2, 1, 0.5) + 0.1,mgp=c(1.75, 0.7, 0))
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
	 if(input$type == "2Dmito") {inpt = "theta (VDAC1)"}else{inpt = input$type}
	 dvs = dat[(dat$type==inpt),]
     mins = aggregate(dvs$value,by=list(dvs$ch),FUN=min)
     maxs = aggregate(dvs$value,by=list(dvs$ch),FUN=max)
     diffs = maxs
     diffs$x = diffs$x - mins$x

     #dscl = d()
	 dscl = updateDat(dat, inpt, input$subject, input$hichan, cord)
	 ids = selected$ids
	 mlab = inpt
	 if(length(ids)>3) {
	   dscl = dscl[dscl$cell_id%in%ids,]
	   mlab = paste(mlab,"Selected fibres",sep="\n") 
	 }
	 
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
	   axrng = range(dscl$value[(is.finite(dscl$value))&(dscl$type==inpt)],na.rm=TRUE)
	 }else{
	   axrng = c(-1,1)
	 }
	 #print(head(dscl))
	 #print(inpt)
	 #print(axrng)
     corrgram(widef,order=FALSE,lower.panel=plotIMC.shadecor,upper.panel=plotIMC.pts,abs=TRUE,col.regions = colorRampPalette(c("blue","white","red")),xlim=axrng,ylim=axrng,main=mlab)	
	}
  
  output$IMC_mainplot = 
    renderPlot({
	  mystripchart()
    },width=950,height=850,pointsize=26)
  
  #output$IMC_arrayplot <- renderPlot({
  #  myarrayplot()
  #},width=850,height=850,pointsize=20)
  
	output$NDUFA13 <- renderPlot({
	op = par(mar=c(4, 2.75, 1, 0.5) + 0.1,mgp=c(1.75, 0.7, 0))
	    arrayplot(d(),ctrld(),ratdat(),cord=cord,ch="NDUFA13",ids=selected$ids,reg_diff = cond(),mitochan="VDAC1",hichan=input$hichan,showControls=input$showControls,axrngCheck = input$axrngCheck, axrng=input$axrng, chlab=paste("A (",chlabs["NDUFA13"],")",sep=""),logify = TRUE)
	par(op)
	})
    output$NDUFB8 <- renderPlot({
	op = par(mar=c(4, 2.75, 1, 0.5) + 0.1,mgp=c(1.75, 0.7, 0))
	    arrayplot(d(),ctrld(),ratdat(),cord=cord,ch="NDUFB8",ids=selected$ids,reg_diff = cond(),mitochan="VDAC1",hichan=input$hichan,showControls=input$showControls,axrngCheck = input$axrngCheck, axrng=input$axrng, chlab=paste("B (",chlabs["NDUFB8"],")",sep=""),logify = TRUE)
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
    filename =  paste("figures","pdf",sep="."),
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
	  if(input$type=="2Dmito"){
	  pdf(file, pointsize=16,width =  8.27, height =  8.27*4/3)
        myarrayplot()
	  }else{
	  pdf(file, pointsize=16,width = 8.27, height = 8.27)
	    mystripchart()
	  }
      dev.off()
    }
  )	

    output$download_png <- downloadHandler(
	
    filename =  paste("figures","png",sep="."),
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
	  if(input$type=="2Dmito"){
       png(file, width=1750, height=1750*4/3, pointsize=51)
        myarrayplot()
	  }else{
       png(file, width=1750, height=1750, pointsize=51)
	    mystripchart()
	  }
      dev.off()  # turn the device off
    }
  )	


}
enableBookmarking("url")
shinyApp(ui,server)
