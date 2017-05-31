library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(GGally)
library(edgeR)
library(plotly)
library(htmlwidgets)
library(shinyBS)
library(shiny)

ui <- shinyUI(pageWithSidebar(
  headerPanel("Click the button"),
  sidebarPanel(
    uiOutput("selInput"), ########## choose all pairs
    uiOutput("slider"),
    sliderInput("threshP", "P-value:", min = 0, max = 1, value=0.1, step=0.05),
    uiOutput("uiExample") #, ########## choose action GoButton
    #uiOutput("testPair")
    #actionButton("goButton", "Go!")
  ),
  mainPanel(
    plotlyOutput("plot1"),
    #verbatimTextOutput("click"),
    #verbatimTextOutput("selectedValues"),
    plotlyOutput("boxPlot")
  )
))

dds <- readRDS("/Users/lindz/BeeVirusDiet/beeDataDDSRLD.rds")[[1]]
rld <- readRDS("/Users/lindz/BeeVirusDiet/beeDataDDSRLD.rds")[[2]]

myLevels <- unique(sapply(colnames(assay(rld)), function(x) unlist(strsplit(x,"[.]"))[1]))
myPairs <- list()

# Runs exact test on all pairs of groups and saves in list
k=1
for (i in 1:(length(myLevels)-1)){
  for (j in (i+1):(length(myLevels))){
    myPairs[[k]] <- paste(myLevels[i], " and ", myLevels[j])
    k=k+1
  }
}

dat <- readRDS("/Users/lindz/BeeVirusDiet/beeVolcanoData.rds")

nCol = ncol(dat)
datFCP = dat[,(nCol-2*length(myLevels)+1):nCol]
# x-axis FC, y-axis pval

xMax = max(unlist(lapply(seq(1,ncol(datFCP),by=2), function(x) max(datFCP[,x], na.rm=TRUE))))






xMin = min(datFCP[,seq(1,ncol(datFCP),by=2)])
yMax = max(datFCP[,seq(2,ncol(datFCP),by=2)])
yMin = min(datFCP[,seq(2,ncol(datFCP),by=2)])
fcMax = ceiling(max(exp(xMax), 1/exp(xMin)))

boxDat <- dat[, 1:(ncol(dat)-2*length(myLevels))] %>% gather(key, val, -c(ID))
BP <- ggplot(boxDat, aes(x = key, y = val)) + geom_boxplot()
ggBP <- ggplotly(BP)

server <- shinyServer(function(input, output) {
  
  #make dynamic slider
  output$slider <- renderUI({
    sliderInput("threshFC", "Fold change:", min=0, max=fcMax, value=ceiling((fcMax)/3), step=0.5)
  })
  
  output$selInput <- renderUI({
    selectInput("selPair", "Pairs:", myPairs)
  })
  
  pairNum <- reactive(as.numeric(which(myPairs==input$selPair)))
  col1 <- reactive(colnames(dat)[nCol-2*length(myLevels)+2*pairNum()])
  col2 <- reactive(colnames(dat)[nCol-2*length(myLevels)+2*pairNum()-1])
  
  #output$testPair <- renderPrint({str(col1())})
  
  # datInput only validated once the go button is clicked
  datInput <- eventReactive(input$goButton, {
    dat[ which(dat[col1()] > -1* log10(input$threshP) & exp(abs(dat[col2()])) > input$threshFC), ]
  })
  
  output$uiExample <- renderUI({
    tags$span(
      tipify(actionButton("goButton", "Go!"), "A Pointless Button", "This button is pointless!")
    )
  })
  
  threshP <- reactive(input$threshP)
  threshFC <- reactive(input$threshFC)
  
  df <- data.frame()
  p <- ggplot(df) + geom_point() + xlim(xMin, xMax) + ylim(yMin, yMax)
  gp <- ggplotly(p)
  
  output$plot1 <- renderPlotly({gp %>% onRender("
    function(el, x, data) {
    
    var dat = data.dat
    var selFC = [];
    var selP = [];
    var sselID = [];
    dat.forEach(function(row){
    rowP = row[data.col1]
    rowFC = row[data.col2]
    rowID = row['ID']
    selFC.push(rowFC);
    selP.push(rowP);
    sselID.push(rowID);
    });
    
    var Traces = [];
    var tracePoints = {
    x: selFC,
    y: selP,
    text: sselID,
    mode: 'markers',
    marker: {
    color: 'black',
    size: 2
    },
    };
    Traces.push(tracePoints);
    Plotly.addTraces(el.id, Traces);
    
    el.on('plotly_selected', function(e) {
    numSel = e.points.length
    Points = e.points
    selID = []
    for (a=0; a<numSel; a++){
    PN = Points[a].pointNumber
    selRow = sselID[PN]
    selID.push(selRow)
    }
    Shiny.onInputChange('selID', selID);
    })
    }", data = list(dat = datInput(), col1 = isolate(col1()), col2 = isolate(col2())))})

  selID <- reactive(input$selID)
  pcpDat <- reactive(dat[which(dat$ID %in% selID()), 1:(ncol(dat)-2*length(myLevels))])
  #output$selectedValues <- renderPrint({str(pcpDat())})
  colNms <- colnames(dat[, 2:(ncol(dat)-2*length(myLevels))])
  nVar <- length(2:(ncol(dat)-2*length(myLevels)))
  
  output$boxPlot <- renderPlotly({
    ggBP %>% onRender("
      function(el, x, data) {
      
      var Traces = [];
      
      var dLength = data.pcpDat.length
      var vLength = data.nVar
      var cNames = data.colNms
      
      for (a=0; a<dLength; a++){
      xArr = [];
      yArr = [];
      for (b=0; b<vLength; b++){
      xArr.push(b+1)
      yArr.push(data.pcpDat[a][cNames[b]]);
      }
      
      var traceHiLine = {
      x: xArr,
      y: yArr,
      mode: 'lines',
      line: {
      color: 'orange',
      width: 1
      },
      opacity: 0.9,
      }
      Traces.push(traceHiLine);
      }
      Plotly.addTraces(el.id, Traces);
      
      }", data = list(pcpDat = pcpDat(), nVar = nVar, colNms = colNms))})

})

shinyApp(ui, server)



