library(shiny)
library(ggplot2)
source('~/Code/NetworkModeling/R/DB/DBFunctions.R')

makeHist = function(dataset){
  con = myConnect()
  df = data.frame()
  for(d in dataset){
    res = dbGetQuery(con, sprintf("select log2fc from %s", dataset))  
    df = rbind(df, data.frame(type=d, log2fc=res$log2fc))
  }
  dbDisconnect(con)
  return(df)
}


# Define server logic required to draw a histogram
shinyServer(
  
  function(input, output) {
  
  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot
  
  output$selectedDatasets <- renderText({ 
    sprintf("You have selected datasets :\n %s ", paste(input$checkGroup, collapse=', '))
  })
  
  output$selectedHist <- renderPlot({
    df = makeHist(input$checkGroup)
    p = ggplot(df, aes(x=log2fc, fill=type))
    p + geom_histogram(alpha=.5)
  })
})