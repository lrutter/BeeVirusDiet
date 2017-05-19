library(plotly)
library(htmlwidgets)
library(shiny)
library(GGally)

ui <- shinyUI(fluidPage(
  sliderInput("ci", "Prediction interval:", min = 0, max = 0.99, value=0.95, step=0.01),
  plotlyOutput("myPlot")
))

server <- shinyServer(function(input, output) {

  ci <- reactive(input$ci)

  beeCounts <-read.delim(file="../AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)
  colnames(beeCounts) <- c("NC.1", "NC.2", "NR.1", "VR.1", "NS.1", "VP.1", "NS.2", "VR.2", "NP.1", "VP.2", "VC.1", "NP.2", "VP.3", "NP.3", "VS.1", "VS.2", "VC.2", "NC.3", "VP.4", "NC.4", "NR.2", "VC.3", "VC.4", "NP.4", "VR.3", "NC.5", "VS.3", "NP.5", "VC.5", "VS.4", "NS.3", "VS.5", "VP.5", "NR.3", "NR.4", "VC.6", "NS.4", "NC.6", "NP.6", "VR.4", "NR.5", "NR.6", "NS.5", "VP.6", "NS.6", "VR.5", "VR.6", "VS.6")
  countdata <- beeCounts[ , order(names(beeCounts))]
  countdata <- as.matrix(countdata)
  coldata = data.frame(row.names = colnames(countdata), virus = unlist(lapply(colnames(countdata), function (x) substring(unlist(strsplit(x, "[.]"))[1],1,1))), diet = unlist(lapply(colnames(countdata), function (x) substring(unlist(strsplit(x, "[.]"))[1],2,2))), treatment = unlist(lapply(colnames(countdata), function (x) unlist(strsplit(x, "[.]"))[1])))
  dds = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ treatment)
  dds <- DESeq(dds)
  rld <- rlog(dds)
  dat <- as.data.frame(assay(rld))
  setDT(dat, keep.rownames = TRUE)[]
  colnames(dat)[1] <- "ID"
  dat$ID <- as.character(dat$ID)
  colNames <- colnames(dat)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  dat <- as.data.frame(dat)

  nCol = ncol(dat)

  conf=seq(0,0.99,.01)
  st<- qt(1-(1-conf)/2,(nrow(dat)-2))

  b0 = c()
  b1 = c()
  sse = c()
  for (i in 2:(nCol-1)){
    j = nCol
    while (j >i){
      datXY <- as.data.frame(cbind(x = dat[,i], y = dat[,j]))
      datLm <- lm(y~x,data=datXY)
      b0 <- c(b0, coef(datLm)[1])
      b1 <- c(b1, coef(datLm)[2])
      sse <- c(sse, summary(datLm)[[6]])
      j = j-1
    }
  }
  b0<-as.vector(b0)
  b1<-as.vector(b1)
  sse<-as.vector(sse)

  minVal = 0
  maxVal = max(dat[,-1])

  my_fn <- function(data, mapping, ...){
    x = data[,c(as.character(mapping$x))]
    y = data[,c(as.character(mapping$y))]
    p <- ggplot(data = dat, aes(x=x, y=y)) + coord_cartesian(xlim = c(minVal, maxVal), ylim = c(minVal, maxVal))
    p
  }

  p <- ggpairs(dat[,-1], lower = list(continuous = my_fn))

  ggPS <- ggplotly(p)

  myLength <- length(ggPS[["x"]][["data"]])
  for (i in 1:myLength){
    item =ggPS[["x"]][["data"]][[i]]$text[1]
    if (!is.null(item))
      if (!startsWith(item, "co")){
        ggPS[["x"]][["data"]][[i]]$hoverinfo <- "none"
      }
  }



  output$myPlot <- renderPlotly(ggPS %>%
  onRender("
   function(el, x, data) {
  myDat = data.dat

function range(start, stop, step){
var a=[start], b=start;
while(b<stop){b+=step;a.push(b)}
return a;
}

   len = Math.sqrt(document.getElementsByClassName('cartesianlayer')[0].childNodes.length);
   AxisNames = [];
   for (i = 1; i < (len+1); i++) {
   AxisNames.push(document.getElementsByClassName('infolayer')[0].childNodes[i].textContent);
   }

   stIndex = Math.round((1-0)/.01*data.ci)
   st = data.st[stIndex]

   var Traces = [];
   var i=0;
   var j=0;
   var k=1;
   var SubPoints = [];
   while ((i*len+k)<=Math.pow((len-1),2)) {
   while ((i+k)<len){
     var x = [];
     var y = [];
     var xTotal = 0;
     var ssx = 0;
     n = data.dat.length; // 100
     for (a=0; a<n; a++){
     xa = data.dat[a][AxisNames[i]]
     x.push(xa)
     y.push(data.dat[a][AxisNames[(len-k)]])
     xTotal+=xa
     }

     var xm = xTotal/n
     for (a=0; a<n; a++){
     ssx+=Math.pow((data.dat[a][AxisNames[i]] - xm),2)
     }

var minX = -1
var maxX = 2*Math.max.apply(null,x)
     //var minX = Math.min.apply(null,x)
     //var maxX = Math.max.apply(null,x)

     var inc = (maxX-minX)/100
     var xv = [];
     var yv = [];
     var se = [];
     var ci = [];
     var uyv = [];
     var lyv = [];
     var a = minX
     while (a < maxX){
     xv.push(a);
     yva = data.b0[j]+data.b1[j]*a;
     // just changed this to have 1+1/n instead of just 1/n
     sea = data.sse[j] * Math.sqrt(1+1/n+Math.pow((a-xm),2)/ssx);
     yv.push(yva);
     se.push(sea);
     ci.push(st*sea);
     uyv.push(yva+st*sea);
     lyv.push(yva-st*sea);
     a+=inc;
     }

     var lwr = [];
     var upr = [];
     var ypred = [];
     var ssea = [];
     var outCI = [];
     var xPoints = [];
     var yPoints = [];
     var keepIndex = []
     for (a=0; a<n; a++){
     xa = data.dat[a][AxisNames[i]]
     // just changed this to have 1+1/n instead of just 1/n
     ssea.push(data.sse[j] * Math.sqrt(1+1/n+Math.pow((xa-xm),2)/ssx))
     ypred.push(data.b0[j]+data.b1[j]*xa)
     lwr.push(ypred[a] - ssea[a]*st)
     upr.push(ypred[a] + ssea[a]*st)
     if (!(y[a]>lwr[a] & y[a]<upr[a])){
     xPoints.push(xa)
     yPoints.push(data.dat[a][AxisNames[(len-k)]])
     keepIndex.push(a)
     }
     }
   SubPoints.push(keepIndex);

     var tracePoints = {
     x: xPoints,
     y: yPoints,
     mode: 'markers',
     marker: {
     color: 'black',
     size: 4
     },
     xaxis: 'x' + (i+1),
     yaxis: 'y' + (i*len+k),
     hoverinfo: 'none'
     };
     var hiLine = {
     x: xv,
     y: uyv,
     mode: 'lines',
     line: {
     color: 'gray',
     width: 1
     },
     xaxis: 'x' + (i+1),
     yaxis: 'y' + (i*len+k),
     opacity: 0.25,
     hoverinfo: 'none'
     };
     var lowLine = {
     x: xv,
     y: lyv,
     mode: 'lines',
     fill: 'tonexty',
     line: {
     color: 'gray',
     width: 1
     },
     xaxis: 'x' + (i+1),
     yaxis: 'y' + (i*len+k),
     opacity: 0.25,
     hoverinfo: 'none'
     };
     Traces.push(tracePoints);
     Traces.push(hiLine);
     Traces.push(lowLine);
     j++;
     k++;
   }
   i++;
   k=1;
   }
   Plotly.addTraces(el.id, Traces);



var idRows = []
for (a=0; a<data.dat.length; a++){
idRows.push(data.dat[a]['ID'])
}

var nseltrace = 0;
el.on('plotly_selected', function(e) {

numSel = e.points.length
cN = e.points[0].curveNumber;

var pointNumbers = [];
for (a=0; a<numSel; a++){
pointNumbers.push(e.points[a].pointNumber)
}

// Determine which subplot was selected
subPlot = (cN - Math.pow(len,2))/3+1

var selData = []
for (a=0; a<pointNumbers.length; a++){
var selDot = SubPoints[subPlot-1][pointNumbers[a]]
selData.push(myDat[selDot])
}


if (nseltrace>0){
Plotly.deleteTraces(el.id,range((-1*len*(len-1)/2),-1,1))
}

var Traces = [];
var i=0;
var k=1;
while ((i*len+k)<=Math.pow((len-1),2)) {
var xArr = [];
for (a=0; a<selData.length; a++){
xArr.push(selData[a][AxisNames[i]])
}
while ((i+k)<len){
  var yArr = [];
  for (a=0; a<selData.length; a++){
  yArr.push(selData[a][AxisNames[(len-k)]])
  }
  var trace = {
  x: xArr,
  y: yArr,
  mode: 'markers',
  marker: {
  color: 'red',
  size: 4
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

nseltrace = nseltrace+1
Plotly.addTraces(el.id, Traces);

})
   }
   ", data = list(dat=dat, b0=b0, b1=b1, sse=sse, st=st, ci=ci())))})

shinyApp(ui, server)
