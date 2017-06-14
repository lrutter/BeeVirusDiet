# This script creates hexbin background based on only those samples and genes. In other words, uses (50,000 * 12) instead of (50,000 * 48)

library(plotly)
library(GGally)
library(hexbin)
library(htmlwidgets)
library(tidyr)
library(shiny)
library(edgeR)
library(EDASeq)
library(dplyr)
library(data.table)
library(ggplot2)
library(DESeq2)

ui <- shinyUI(fluidPage(
  #uiOutput("selInput"),
  #actionButton("goButton", "Go!"),
  #verbatimTextOutput("selectedValues"),
  plotlyOutput("scatMatPlot")
))

server <- shinyServer(function(input, output) {
  
  dds <- readRDS("/Users/lindz/BeeVirusDiet/tblshoot-VSNP/Method1/beeDataDDSRLD.rds")[[1]]
  rld <- readRDS("/Users/lindz/BeeVirusDiet/tblshoot-VSNP/Method1/beeDataDDSRLD.rds")[[2]]
  
output$scatMatPlot <- renderPlotly({  
  
  # Change group1 and group2 as needed
  group1 ="VP"
  group2 ="VR"
  
  bindataSel <- as.data.frame(assay(rld))[, sampleIndex]
  setDT(bindataSel, keep.rownames = TRUE)[]
  colnames(bindataSel)[1] <- "ID"
  bindataSel$ID <- as.character(bindataSel$ID)
  bindataSel <- as.data.frame(bindataSel) #15314*13
  
  sampleIndex <- which(sapply(colnames(assay(rld)), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(group1, group2))
  sampleIndex1 <- which(sapply(colnames(bindataSel), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(group1))
  sampleIndex2 <- which(sapply(colnames(bindataSel), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(group2))
  
  res <- results(dds, contrast=c("treatment",group1,group2))
  degIndex <- which(res@listData$padj<0.05) 
  resSort <- res[ order(res[,6]), rm.NA=TRUE]
  g <- rownames(resSort)[1]
  currGene <- bindataSel[which(bindataSel$ID==g),]
  currGene1 <- unname(unlist(currGene[,sampleIndex1]))
  currGene2 <- unname(unlist(currGene[,sampleIndex2]))
  len <- length(currGene1)
  
  c1 <- c()
  c2 <- c()
  
  k=1
  for (i in 1:len){
    for (j in 1:len){
      c1[k] <- currGene1[i]
      c2[k] <- currGene2[j]
      k <- k+1
    }
  }
  
  orangeDots <- data.frame(x=c1, y=c2)
  
  minVal = min(bindataSel[,-1])
  maxVal = max(bindataSel[,-1])
  maxRange = c(minVal, maxVal)
  xbins=10
  buffer = (maxRange[2]-maxRange[1])/xbins/2 # Because usually shows at least half of hex?
  
  x = unlist(bindataSel[,(sampleIndex1)])
  y = unlist(bindataSel[,(sampleIndex2)])
  h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
  hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
  attr(hexdf, "cID") <- h@cID
  p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(maxRange[1]-1*buffer, maxRange[2]+buffer), ylim = c(maxRange[1]-1*buffer, maxRange[2]+buffer))
  
  #p + geom_point(aes(x=x, y=y), data=orangeDots, inherit.aes = FALSE, color="orange")
  
  c1
  c2
  
    ggPS <- ggplotly(p)

    myLength <- length(ggPS[["x"]][["data"]])
    for (i in 1:myLength){
      hexHover = ggPS[["x"]][["data"]][[i]]$text
      if (!is.null(hexHover) && grepl("hexID", hexHover)){
        ggPS[["x"]][["data"]][[i]]$text <- strsplit(hexHover, "<")[[1]][1]
        ggPS[["x"]][["data"]][[i]]$t2 <- hexHover
        ggPS[["x"]][["data"]][[i]]$hoverinfo <- "text"
      }
    }
    
    ggPS %>% onRender("
      function(el, x, data) {

      // Save selected row IDs for PCP
      //Shiny.onInputChange('selID', selID);
      console.log('here!')
console.log(data.c1)
console.log(data.c2)
      var Traces = [];
      var trace = {
      x: data.c1,
      y: data.c2,
      mode: 'markers',
      marker: {
      color: 'orange',
      size: 6
      },
      hoverinfo: 'none'
      };
      Traces.push(trace);
      Plotly.addTraces(el.id, Traces);
      }
      ", data = list(c1=c1, c2=c2))
  })
})

shinyApp(ui, server)