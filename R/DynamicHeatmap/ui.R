library(shiny)
source('~/Code/NetworkModeling/R/DB/DBFunctions.R')

getDataSets = function(){
  con <- myConnect()
  res = dbGetQuery(con, "select distinct omics_type from experiments")
  res = as.list(t(res))
  dbDisconnect(con)
  return(res)
}

datasets = getDataSets()
datasets = list('metabolomics_mouse','genomics_rnaseq_mouse','proteomics_ph_sites_mouse','proteomics_ub_sites_mouse')

shinyUI(fluidPage(
  titlePanel("Fold change distributions"),
  
  fluidRow(

    column(3, 
           checkboxGroupInput("checkGroup", 
                              label = h5("H1N1 datasets"), 
                              datasets,
                              selected = datasets[1:length(datasets)]),
           submitButton("Submit"))
  ),


  mainPanel(
    textOutput("selectedDatasets"),
    plotOutput("selectedHist")
  )

  
))