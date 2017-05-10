library(plotly)
library(GGally)
library(hexbin)
library(htmlwidgets)
library(tidyr)
library(shiny)
library(dplyr)
library(data.table)
library(ggplot2)

ui <- shinyUI(fluidPage(
  plotlyOutput("scatMatPlot", height = 700),
  #verbatimTextOutput("selectedValues"),
  plotlyOutput("boxPlot")
))

server <- shinyServer(function(input, output, session) {
  
  set.seed(1)
  bindata <- data.frame(ID = paste0("ID",1:100), A.1=abs(rnorm(100)), A.2=abs(rnorm(100)))
  bindata$ID <- as.character(bindata$ID)
  colNames <- colnames(bindata)
  
  output$scatMatPlot <- renderPlotly({

  maxVal = max(abs(bindata[,-1]))
  maxRange = c(0, maxVal)
  xbins=10
  buffer = maxRange[2]/xbins
  
  my_fn <- function(data, mapping, ...){
    x = data[,c(as.character(mapping$x))]
    y = data[,c(as.character(mapping$y))]
    h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
    hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
    attr(hexdf, "cID") <- h@cID
    p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-1*buffer, maxRange[2]+buffer), ylim = c(-1*buffer, maxRange[2]+buffer))
    p
  }
  
  p <- ggpairs(bindata[,-1], lower = list(continuous = my_fn))
  pS <- p
  
  ggPS <- ggplotly(pS)
  
  myLength <- length(ggPS[["x"]][["data"]])
  for (i in 1:myLength){
    item =ggPS[["x"]][["data"]][[i]]$text[1]
    if (!is.null(item)){
      if (!startsWith(item, "co")){
        ggPS[["x"]][["data"]][[i]]$hoverinfo <- "none"
      }}
    hexHover = ggPS[["x"]][["data"]][[i]]$text
    if (!is.null(hexHover) && grepl("hexID", hexHover)){
      ggPS[["x"]][["data"]][[i]]$text <- strsplit(hexHover, "<")[[1]][1]
      ggPS[["x"]][["data"]][[i]]$t2 <- hexHover
      ggPS[["x"]][["data"]][[i]]$hoverinfo <- "text"
    }
  }
  
  for(i in 2:(p$nrow)) {
    for(j in 1:(p$nrow-1)) {
      bindata[[paste(i,j,sep="-")]] <- attr(pS[i,j]$data, "cID")
    }
  }

  ggPS %>% onRender("
    function(el, x, data) {
    noPoint = x.data.length;
    
    el.on('plotly_click', function(e) {
    
    if (x.data.length > noPoint){
    //Plotly.deleteTraces(el.id, range(noPoint, (noPoint+(len*(len-1)/2-1)), 1));
    }
    
    myX = 2
    myY = 1
    cN = e.points[0].curveNumber
    split1 = (x.data[cN].text).split(' ')
    hexID = (x.data[cN].t2).split(':')[2]
    var selRows = [];
    data.forEach(function(row){
    if(row[myX+'-'+myY]==hexID) selRows.push(row);
    });
    selID = []
    for (a=0; a<selRows.length; a++){
    selID.push(selRows[a]['ID'])
    }
    // Save selected row IDs for PCP
    Shiny.onInputChange('selID', selID);
    
    AxisNames = ['A.1', 'A.2']

    var Traces = [];
    var xArr = [];
    for (a=0; a<selRows.length; a++){
    xArr.push(selRows[a]['A.1'])
    }
    var yArr = [];
    for (a=0; a<selRows.length; a++){
    yArr.push(selRows[a]['A.2'])
    }
    var trace = {
    x: xArr,
    y: yArr,
    mode: 'markers',
    marker: {
    color: 'orange',
    size: 6
    },
    hoverinfo: 'none'
    };
    Traces.push(trace);
    Plotly.addTraces(el.id, Traces);
    })}
    ", data = bindata)
  })
  
  selID <- reactive(input$selID)
  
  output$boxPlot <- renderPlotly({
    
    bindata[2:ncol(bindata)] = log(bindata[2:ncol(bindata)]/colMeans(bindata[2:ncol(bindata)]))
    nVar = ncol(bindata)
    pcpDat <- reactive(bindata[which(bindata$ID %in% selID()), c(1:nVar)])
    
    observe({
      session$sendCustomMessage(type = "lines", pcpDat())
    })
    
    colNms <- colnames(bindata[, c(2:nVar)])
    
    boxDat <- bindata[, c(1:nVar)] %>% gather(key, val, -c(ID))
    colnames(boxDat) <- c("ID", "Sample", "Count")
    
    output$selectedValues <- renderPrint({str(boxDat)})
    
    BP <- ggplot(boxDat, aes(x = Sample, y = Count)) + geom_boxplot()
    ggBP <- ggplotly(BP)
    
    ggBP %>% onRender("
      function(el, x, data) {
      console.log(x.data)
      var Traces =[];
      Shiny.addCustomMessageHandler('lines',
      function(drawLines) {
      console.log(x.data.length)
      if (x.data.length > 1){
        //Plotly.deleteTraces(el.id,1); // not working
        //Plotly.deleteTraces(el.id, x.data.length-1) // not working
      }
      
      var vLength = Object.keys(drawLines).length
      var AxisNamesKeys = Object.keys(drawLines);
      var cNames = [];
      for (i = 1; i < vLength; i++) {
      cNames.push(AxisNamesKeys[i])
      }
      var first = cNames[0]
      dLength = drawLines[first].length      

      for (a=0; a<dLength; a++){
      xArr = [];
      yArr = [];
      for (b=0; b<(vLength-1); b++){
      xArr.push(b+1)
      yArr.push(drawLines[cNames[b]][a]);
      }
      
      var traceLine = {
      x: xArr,
      y: yArr,
      mode: 'lines',
      line: {
      color: 'orange',
      width: 1
      },
      opacity: 0.9,
      }
      Traces.push(traceLine);
      }
      Plotly.addTraces(el.id, Traces);
    })
  }" )})
})

shinyApp(ui, server)