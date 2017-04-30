ui <- shinyUI(fluidPage(
  actionButton("addButton", "Add points"),
  plotlyOutput("myPlot", height = 700)
))

server <- shinyServer(function(input, output) {
  
  output$myPlot <- renderPlotly({
    p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank()
    pP <- ggplotly(p)
    
    pP %>% onRender("
    function(el, x, data) {
$('#addButton').on('click',function() {
      var Traces = [];
      var trace = {
      x: data.x,
      y: data.y,
      mode: 'markers',
      marker: {
      color: 'green',
      size: 6
      },
      hoverinfo: 'none'
      };
      Traces.push(trace);
      Plotly.addTraces(el.id, Traces);
})

    }", data = list(x = mtcars$wt, y = mtcars$mpg))
  })      
})

shinyApp(ui, server)