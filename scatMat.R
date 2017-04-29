library(plotly)
library(data.table)
library(GGally)
library(hexbin)
library(htmlwidgets)
library(tidyr)
library(shiny)
library(edgeR)
library(EDASeq)
library(dplyr)
library(data.table)

ui <- shinyUI(fluidPage(
  sidebarLayout(
    sidebarPanel(
      uiOutput("selInput"),
      uiOutput("go")
    ),
    mainPanel(
      verbatimTextOutput("test"),
      verbatimTextOutput("test2"),
      plotlyOutput("scatMatPlot", height = 700)
      #plotlyOutput("boxPlot"),
      #verbatimTextOutput("selectedValues")
    )))
)

server <- shinyServer(function(input, output) {

  set.seed(1)
  bindata <- data.frame(ID = paste0("ID",1:100), A.1=rnorm(100), A.2=rnorm(100), A.3=rnorm(100), B.1=rnorm(100), B.2=rnorm(100), C.1=rnorm(100), C.2=rnorm(100))
  bindata$ID <- as.character(bindata$ID)
  colNames <- colnames(bindata)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
    
  output$selInput <- renderUI({
    checkboxGroupInput("selPair", "Choose groups:", myPairs)
  })
  
  output$test <- renderPrint({
    str(input$selPair)
  })
  
  output$go <- renderUI({
    tags$span(
      tipify(actionButton("goButton", "Go!"), "Message", "Sub-Message")
    )
  })

  #datInput only validated once the go button is clicked
  datInput <- eventReactive(input$goButton, {
    cbind(ID=bindata$ID, bindata[,which(colGroups %in% input$selPair)])
  })
  
  output$test2 <- renderPrint({
    str(datInput())
  })
  
  output$scatMatPlot <- renderPlotly({

    # will wait to render until datInput is validated
    bindataSel <- datInput()
    maxVal = max(abs(bindataSel[,-1]))
    maxRange = c(-1*maxVal, maxVal)
    #maxRange = c(0, maxVal)

    my_fn <- function(data, mapping, ...){
      x = data[,c(as.character(mapping$x))]
      y = data[,c(as.character(mapping$y))]
      h <- hexbin(x=x, y=y, xbins=10, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
      hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
      attr(hexdf, "cID") <- h@cID
      p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-0.5, maxRange[2]+0.5), ylim = c(-0.5, maxRange[2]+0.5))
      p
    }

    p <- ggpairs(bindataSel[,-1], lower = list(continuous = my_fn))
    pS <- p

    ggPS <- ggplotly(pS)
    #ggPS <- ggplotly(p)

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
        bindataSel[[paste(i,j,sep="-")]] <- attr(pS[i,j]$data, "cID")
      }
    }


    ggPS %>% onRender("
    function(el, x, data) {

    function range(start, stop, step){
    var a=[start], b=start;
    while(b<stop){b+=step;a.push(b)}
    return a;
    };

    len = Math.sqrt(document.getElementsByClassName('cartesianlayer')[0].childNodes.length);
    AxisNames = [];
    for (i = 1; i < (len+1); i++) {
    AxisNames.push(document.getElementsByClassName('infolayer')[0].childNodes[i].textContent);
    }
    noPoint = x.data.length;

    el.on('plotly_click', function(e) {
    console.log(e)

    if (x.data.length > noPoint){
    Plotly.deleteTraces(el.id, range(noPoint, (noPoint+(len*(len-1)/2-1)), 1));
    }

    xVar = (e.points[0].xaxis._id).replace(/[^0-9]/g,'')
    if (xVar.length == 0) xVar = 1
    yVar = (e.points[0].yaxis._id).replace(/[^0-9]/g,'')
    if (yVar.length == 0) yVar = 1
    myX = len + 1 - (yVar - len * (xVar - 1))
    console.log(['myX', myX])
    myY = xVar
    console.log(['myY', myY])
    cN = e.points[0].curveNumber
    //console.log(['cN', cN])
    split1 = (x.data[cN].text).split(' ')
    console.log(['split1',split1])
    console.log(x.data[cN].t2)
    hexID = (x.data[cN].t2).split(':')[2]
    console.log(['hexID', hexID])
    counts = split1[1].split('<')[0]
    var selRows = [];
    data.forEach(function(row){
    if(row[myX+'-'+myY]==hexID) selRows.push(row);
    });
    console.log(selRows)
    selID = []
    for (a=0; a<selRows.length; a++){
    selID.push(selRows[a]['ID'])
    }
    console.log(selID)
    // Save selected row IDs for PCP
    Shiny.onInputChange('selID', selID);

    var Traces = [];
    var i=0;
    var k=1;
    while ((i*len+k)<=Math.pow((len-1),2)) {
    var xArr = [];
    for (a=0; a<selRows.length; a++){
    xArr.push(selRows[a][AxisNames[i]])
    }
    while ((i+k)<len){
    var yArr = [];
    for (a=0; a<selRows.length; a++){
    yArr.push(selRows[a][AxisNames[(len-k)]])
    }
    var trace = {
    x: xArr,
    y: yArr,
    mode: 'markers',
    marker: {
    color: 'orange',
    size: 6
    },
    xaxis: 'x' + (i+1),
    yaxis: 'y' + (i*len+k),
    hoverinfo: 'none'
    };
    Traces.push(trace);
    k++;
    }
    i++;
    k=1;
    }
    Plotly.addTraces(el.id, Traces);
    })}
    ", data = bindataSel)
  })

  #   selID <- reactive(input$selID)
  # 
  # pcpDat <- reactive(bindata[which(bindata$ID %in% selID()), c(1:(p$nrow+1))])
  # output$selectedValues <- renderPrint({str(pcpDat())})
  # colNms <- colnames(bindata[, c(2:(p$nrow+1))])
  # 
  # boxDat <- bindata[, c(1:(p$nrow+1))] %>% gather(key, val, -c(ID))
  # BP <- ggplot(boxDat, aes(x = key, y = val)) + geom_boxplot()
  # ggBP <- ggplotly(BP)
  
  # output$boxPlot <- renderPlotly({
  #   ggBP %>% onRender("
  #   function(el, x, data) {
  #   
  #   var Traces = [];
  #   
  #   var dLength = data.pcpDat.length
  #   var vLength = data.nVar
  #   var cNames = data.colNms
  #   
  #   for (a=0; a<dLength; a++){
  #   xArr = [];
  #   yArr = [];
  #   for (b=0; b<vLength; b++){
  #   xArr.push(b+1)
  #   yArr.push(data.pcpDat[a][cNames[b]]);
  #   }
  #   
  #   var traceHiLine = {
  #   x: xArr,
  #   y: yArr,
  #   mode: 'lines',
  #   line: {
  #   color: 'orange',
  #   width: 1
  #   },
  #   opacity: 0.9,
  #   }
  #   Traces.push(traceHiLine);
  #   }
  #   Plotly.addTraces(el.id, Traces);
  #   
  #   }", data = list(pcpDat = pcpDat(), nVar = p$nrow, colNms = colNms))})
})

shinyApp(ui, server)


