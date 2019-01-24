library(shiny)
library(data.table)

dat = fread("../dat.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
dat$ch = substring(dat$channel,regexpr("\\_[^\\_]*$", dat$channel)+1,nchar(dat$channel))
dat$type = "Mean intensity"
dat$type[grepl("LOG_",dat$channel)] = "Log mean intensity"
dat$type[grepl("MED_",dat$channel)] = "Median intensity"
dat$type[grepl("R_",dat$channel)] = "Ratio mean intensity (VDAC1)"
dat$type[grepl("R_MED_",dat$channel)] = "Ratio median intensity (VDAC1)"
dat$type[grepl("R_LOG_",dat$channel)] = "Ratio log mean intensity (VDAC1)"
dat$type[grepl("Z_",dat$channel)] = "z-score"
types = unique(dat$type)

subdf = unique(dat[,c("patient_id","subject_group")])
subjs = subdf$patient_id
grps = subdf$subject_group
names(grps) = subjs
labs = paste(subjs," (",grps[subjs],")",sep="")
names(subjs) = labs

cutcords = c(2.5,3.5,4.5,6.5,7.5,8.5)
cordlabs = c("CI","CII","CIII","CIV","CV","OMM","Cell")
cord = c("NDUFB8","GRIM19","SDHA","UqCRC2","COX4","MTCO1","OSCP","VDAC1","Dystrophin","DNA1")
dat = dat[dat$ch%in%cord,]
dat$ch = factor(dat$ch, levels = cord)

sgroups = list(
ctrl = c("C_m1180","C_M1089","C_M1217"),
CI = c("P_M0838","P_M0517"), 
large = c("P_M0966","P_M0284"),
"MT-TL1" = c("P_M0694","P_m0917","P_M963"),
tRNA = c("P_M0913","P_M1513","P_M0207")
)
sgrouplist = unlist(sgroups)

group_order =c("ctrl","CI","large","MT-TL1","tRNA")
dat = dat[order(match(dat$subject_group,group_order),decreasing = FALSE),]
dat$patient_id = factor(dat$patient_id,unique(dat$patient_id))
dat$subject_group = factor(dat$subject_group,unique(dat$subject_group))
subtab = aggregate(patient_id~subject_group,data=dat,function(x) length(unique(x)))
nsubj = subtab$patient_id

nums = 1:length(sgrouplist)
names(nums) = sgrouplist

dat$num = nums[dat$patient_id]
dat$jit = dat$num + runif(length(dat$num),0,0.3)

nmax = 20

ui <- fluidPage(
  titlePanel("Warren et al. (2018) IMC results (by channel)"),
  sidebarLayout(
    sidebarPanel(
	width = 2,
      selectInput("channel", "Channel", cord),
	selectInput("type", "Measure", types, selected="z-score"),
	p("Drag a box on the plot area to select fibres.")
    ),
   mainPanel(
    plotOutput("IMC_stripchart",height=1000,
	click = "plot_click",
      hover = hoverOpts(id = "plot_hover", delayType = "throttle"),
      brush = brushOpts(id = "plot_brush")
	),
        h4("Selected fibres"),
        tableOutput("plot_brushedpoints")
    )
  )
)


server <- function(input, output) {
  vals = reactiveValues(
    pts = dat[dat$channel=="NOCHANNEL",]
  )

  output$plot_brushedpoints = renderTable(vals$pts)
  

  output$IMC_stripchart <- renderPlot({
    
    dt = dat[dat$type==input$type,]
    d = dt[dt$ch==input$channel,]
    vals$pts = brushedPoints(d, input$plot_brush, "jit", "value")
    
    subjgrp = unique(d$ch)
    main = input$channel

    op = par(mar=c(6, 4, 2, 0) + 0.1)
    if(grepl("Ratio",input$type)){logval="y"}else{logval=""}
    plot(d$jit,d$value,xlab="",ylab = paste(input$channel,input$type),main = main,col = rgb(0,0,0,0.05), pch = 16, cex=0.5,las=2,axes=FALSE,log=logval)
    axis(2)
    axis(1,at=1:length(nums),labels=names(nums),las=2)
    abline(v=cumsum(nsubj[1:(length(nsubj)-1)])+0.5,lty=2)
    text(c(0,cumsum(nsubj[1:(length(nsubj)-1)]))+0.25,max(d$value),as.character(subtab$subject_group),pos=4)
    ids = unique(vals$pts$cell_id)
    #if(length(ids)>nmax) ids = sample(ids,nmax)
    hl = d[d$cell_id%in%ids,]
    hl = hl[order(hl$jit),]
    for(i in ids){
	points(hl$jit[hl$cell_id==i],hl$value[hl$cell_id==i],type="b",lwd=2,col=rgb(1,0,0,0.3))
    }
    par(op)
  },width=1000,height=1000,pointsize=30)
}

shinyApp(ui,server)

