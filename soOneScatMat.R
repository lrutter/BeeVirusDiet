library(plotly)
library(GGally)
library(hexbin)
library(htmlwidgets)
library(tidyr)
library(shiny)
library(dplyr)
library(data.table)
library(ggplot2)
library(tibble)

myPairs <- c("A", "B", "C", "D")

ui <- shinyUI(fluidPage(
  titlePanel("title panel"),
  
  sidebarLayout(position = "left",
    sidebarPanel(
      selectizeInput("selPair", "Pairs:", choices = myPairs, multiple = TRUE, options = list(maxItems = 2)),
      actionButton("goButton", "Go!"),
      width = 3
    ),
    mainPanel(
      verbatimTextOutput("info"),
      plotlyOutput("scatMatPlot")
    )
  )
))

server <- shinyServer(function(input, output, session) {
  
  # Create data and subsets of data based on user selection of pairs
  dat <- data.frame(ID = paste0("ID", 1:10000), A = rnorm(10000), B = rnorm(10000), C = rnorm(10000), D = rnorm(10000))
  pairNum <- reactive(input$selPair)
  group1 <- reactive(pairNum()[1])
  group2 <- reactive(pairNum()[2])
  sampleIndex <- reactive(which(colnames(dat) %in% c(group1(), group2())))
  
  # Create data subset based on two letters user chooses
  datSel <- eventReactive(sampleIndex(), {
    datSel <- dat[, c(1, sampleIndex())]
    datSel$ID <- as.character(datSel$ID)
    datSel <- as.data.frame(datSel)
    datSel
  })
  
  sampleIndex1 <- reactive(which(colnames(datSel()) %in% c(group1())))
  sampleIndex2 <- reactive(which(colnames(datSel()) %in% c(group2())))

  # Create background Plotly graph with hex binning all 100 rows of the two user-selected columns
  ggPS <- eventReactive(datSel(), {
    minVal = min(datSel()[,-1])
    maxVal = max(datSel()[,-1])
    maxRange = c(minVal, maxVal)
    xbins=7
    buffer = (maxRange[2]-maxRange[1])/xbins/2
    x = unlist(datSel()[,(sampleIndex1())])
    y = unlist(datSel()[,(sampleIndex2())])
    h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
    hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
    attr(hexdf, "cID") <- h@cID
    p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(maxRange[1]-1*buffer, maxRange[2]+buffer), ylim = c(maxRange[1]-1*buffer, maxRange[2]+buffer)) + coord_equal(ratio=1) + labs(x = colnames(datSel()[sampleIndex1()]), y = colnames(datSel()[sampleIndex2()]))
    ggPS <- ggplotly(p)
    ggPS})

  # Output hex bin plot created just above
  output$scatMatPlot <- renderPlotly({
    # Each time user pushes Go! button, the next row of the data frame is selected
    datInput <- eventReactive(input$goButton, {
      g <- datSel()$ID[input$goButton]
      
      # Output ID of selected row
      output$info <- renderPrint({
        g
      })
      
      # Get x and y values of seleced row
      currGene <- datSel()[which(datSel()$ID==g),]
      currGene1 <- unname(unlist(currGene[,sampleIndex1()]))
      currGene2 <- unname(unlist(currGene[,sampleIndex2()]))
      c(currGene1, currGene2)
    })
    
    # Send x and y values of selected row into onRender() function
    observe({
      session$sendCustomMessage(type = "points", datInput())
    })

    # Use onRender() function to draw x and y values of seleced row as orange point
    ggPS() %>% onRender("
      function(el, x, data) {

      noPoint = x.data.length;

      Shiny.addCustomMessageHandler('points', function(drawPoints) {
        if (x.data.length > noPoint){
          Plotly.deleteTraces(el.id, x.data.length-1);
        }
  
        var Traces = [];
        var trace = {
          x: drawPoints.slice(0, drawPoints.length/2),
          y: drawPoints.slice(drawPoints.length/2, drawPoints.length),
          mode: 'markers',
          marker: {
            color: 'orange',
            size: 7
          },
          hoverinfo: 'none'
        };
        Traces.push(trace);
        Plotly.addTraces(el.id, Traces);
      });}")
    })
})

shinyApp(ui, server)