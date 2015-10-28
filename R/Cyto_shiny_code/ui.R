
#setwd("~/Box Sync/Doug/projects/FLU_Networking_code/Cyto_shiny_dev_generic_data_mult_2")
#devtools::install_github('rstudio/DT')

#install.packages("shiny", repos=c("http://rstudio.org/_packages", "http://cran.rstudio.com")) 
#update.packages(ask = FALSE)
library(shiny)
library(shinythemes)
library('data.table')
library('DT')
#library(tidyr)
#library('dplyr')
library('reshape2')
#this function takes r dataframes of nodes and edges and outputs the info/data in a list of 
# character strings which is then used by cytoscape.js to create the nework graph. 
source("./Linked_subdirectories/R_to_Cyjs_wrapper_functions/CCN2_pie.R")
# Cytoscape 'wrapper' function, uses output from CCN2_pie.R to create the cytoscape graph.
source("./Linked_subdirectories/R_to_Cyjs_wrapper_functions/cytoscapeJsSimpleNetwork3_strain.R")
# These are R functions used here in this shiny app but generic enough for non-shiny networking applications their uses 
# are found in "Create_Flu_nework#.R"    
source("./Linked_subdirectories/R_functions_2/FUNC_for_example_genomics3.R")
# all data is read in here.
source("Linked_subdirectories/data/model_data_read.R")
# This is the workhorse R script where the node and edge selections are made along with attributes, colors, shapes. Output from here is 
# used by CCN2_pie.R.
source("./Linked_subdirectories/R_functions_2/Create_Flu_network6.R")
#load("./Linked_subdirectories/data/human_rnaseq_data")

# DF1 <<- get_time(genomics_mouse)
# time_points <<- sort(unique(DF1$time))



# These will be needed when subnetworking is done.  
#library(BioNet)
#library(DLBCL)
#library(igraph)
# themes = cerulean, flatly, default, journal, readable, spacelab, united
#(theme = "bootstrap.css",
shinyUI(fluidPage(theme = shinytheme("readable"), 
                  tags$style(type="text/css",
                             ".shiny-output-error { visibility: hidden; }",
                             ".shiny-output-error:before { visibility: hidden; }"
                  ),
                  tags$head(
                    #  Any javascript sourced in here must be placed in the "www" subdirectory in order for it to be found
                    #  remember to add commas if more than one "tag", last one has no comma!
                    tags$script(src="cytoscape.min.js")
                    #  tags$script(src = 'http://cytoscape.github.io/cytoscape.js/api/cytoscape.js-latest/cola.v3.min.js'), 
                  ),
                  tags$head(
                    tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")
                  ),
                  
                  # Application title
                  #    titlePanel(textOutput('title')),
                  
                  sidebarLayout(
                    sidebarPanel(
                      tags$head(
                        #                         tags$style(type="text/css", "select { width: 400px; }"),
                        tags$style(type="text/css", "select { width: 100px; }"),
                        tags$style(type="text/css", "textarea { max-width: 185px; }"),
                        tags$style(type="text/css", ".jslider { max-width: 100px; }"),
                        tags$style(type='text/css', ".well { max-width: 330px; }"),
                        tags$style(type='text/css', ".span4 { max-width: 230px; }")
                      ),
                      width=3, 
                      selectInput("selectSPECIES", label = h6("SPECIES"),
                                  choices = list("HUMAN" = 'human', "MOUSE" = 'mouse'), 
                                  selected = 'mouse'),
                      # This outputs the dynamic UI component
                      uiOutput("ui1"),
                      ###### new stuff
                      uiOutput("ui2"),
                      ###### old stuff
#                       selectInput("selectOMICS", label = h6
#                               ("OMICS-type"),
#                                   choices = list("PROTEOMICS" = 'proteomics', "GENOMICS" = 'genomics',
#                                                  'METABOLOMICS' = 'metabolomics'), 
#                                   selected = 'genomics', multiple = T),
                  
                      sliderInput("q_value", "primary QVAL:", min=.001, max=.1, value=.05),
                      sliderInput("log2FC", "Fold Change:", min=1.5 , max=6.0, value=2.0),
                      div( style="display:inline-block ; color:#FF0000" ,
                           h6("Select virus(es) and timepoints, lighter node color = later time", align = "left")),
                      # The 'div(style="display:inline-block",)' portion simply allows the checkboxes to be displayed side by side, otherwise a new line will start for
                      # each line here.  
                      div( style="display:inline-block ; color:#FF0000" ,checkboxInput("H1N1_s", "H1N1", value = T)),
                      div( style="display:inline-block ; color:#0000FF",checkboxInput("H3N2_s", "H3N2", value = F)),
                      div(style="display:inline-block ; color:#006600", checkboxInput("H5N1_s", "H5N1", value = F)),  
                      div(style="display:inline-block", checkboxInput("SOV", 
                          "Show unselected viruses that also occur on nodes in this network?", value = F)),  
                      if(6 > 0) {
                        div(style="display:inline-block", checkboxInput("T1", textOutput('TP1') , value = F))
                      } ,# }), # close reactive
                      if(6 > 1) {
                        div(style="display:inline-block", checkboxInput("T2", textOutput('TP2'), value = T))
                      },
                      if(6 > 2) {
                        div(style="display:inline-block" , checkboxInput("T3",  textOutput('TP3'), value = F))
                      },
                      if(6 > 3) {
                        div(style="display:inline-block" , checkboxInput("T4", textOutput('TP4'), value = F))
                      },
                      if(6 > 4) {
                        div(style="display:inline-block" , checkboxInput("T5", textOutput('TP5'), value = F))
                      },
                      #                      if(6 > 5) {
                      #                              div(style="display:inline-block", checkboxInput("T6", output.Time_Points()[6], value = F))
                      #                      },
                      #                      if(length(time_points) > 6) {
                      #                              div(style="display:inline-block" , checkboxInput("T7", output.Time_Points()[7], value = F))
                      #                      },
                      #                      if(length(time_points) > 7) {
                      #                              div(style="display:inline-block" , checkboxInput("T8", output.Time_Points()[8], value = F))
                      #                      },    
                      
                      div(style="display:inline-block",checkboxInput("RN_ne", "Remove Nodes without Edges", value = TRUE)),
                      
                      
                      #      submitButton("Update Network"),
                      #            checkboxInput("addLinks", "Add Links on Nodes?", TRUE)
                      div(style="display:inline-block",
                          h5("Uncheck to de-select Edge types: ", align = "left")),
                      # c("literature","metabolic","binary","regulatory","complexes","kinase","signaling")
                      #   color = c("#FF0000", "#000000", "#0000FF", "#008000", "#800000","#8D38C9" ,"#F88017")
                      div(style="display:inline-block", 
                          conditionalPanel(
                            condition = "output.litL",
                            div(style="display:inline-block ;color: #FF0000", checkboxInput("lit", "literature", value = T))
                            # uiOutput('fullUI')
                          )),
                      div(style="display:inline-block",
                          conditionalPanel(
                            condition = "output.metL",
                            div(style= paste0("color:", "#000000"), checkboxInput("met", "metabolic", value = T))
                            #uiOutput('fullUImet')
                          )),
                      #   uiOutput('fullUI'),
                      div(style="display:inline-block",                            
                          conditionalPanel(
                            condition = "output.binL",
                            #       uiOutput('fullUI')
                            div(style="display:inline-block ;color: #0000FF", checkboxInput("bin", "binary", value = T))
                          )),
                      div(style="display:inline-block",
                          conditionalPanel(
                            condition = "output.regL",
                            div(style="display:inline-block ; color:#008000", checkboxInput("reg", "regulatory", value = T))
                          )),
                      div(style="display:inline-block",
                          conditionalPanel(
                            condition = "output.comL",
                            div(style="display:inline-block ;color:#800000", checkboxInput("com", "complexes", value = T))
                          )),
                      div(style="display:inline-block",
                          conditionalPanel(
                            condition = "output.kinL",
                            div(style="display:inline-block ;color:#8D38C9", checkboxInput("kin", "kinase", value = T))
                          )),
                      div(style="display:inline-block",
                          conditionalPanel(
                            condition = "output.sigL",
                            div(style="display:inline-block ;color:#F88017", checkboxInput("sig", "signaling", value = T))
                          )),
                      #                      verbatimTextOutput('x4'),
                      selectInput("layout", "Layout:", 
                                  #                         choices=c("grid", "random", "circle", "breadthfirst", "cose",
                                  #                                   "springy", "cola", "dagre", "force-directed"), 
                                  choices=c("grid", "random", "circle", "cose"), 
                                  
                                  selected="cose", 
                                  multiple=FALSE)
                      
                    ),
                    
                    mainPanel( # div(id="cy"),
                      tabsetPanel(
                        #                  tabPanel("CytoscapeJS_2", div(id="cy"), htmlOutput("cytoscapeJsPlot_2") ),
                        tabPanel("Nodes", DT::dataTableOutput("cytoscapeJsTable_nodes")),
                        #                  tabPanel("Selected Nodes", DT::dataTableOutput("cytoscapeJsTable_nodes2")),
                        tabPanel("CytoscapeJS", div(id="cy"), htmlOutput("cytoscapeJsPlot") ),
                        tabPanel("Nodes_attr", DT::dataTableOutput("cytoscapeJsTable_nodes_attributes")),
                        tabPanel("Edges", DT::dataTableOutput("cytoscapeJsTable_edges")),
                        tabPanel("Help: Select",  includeMarkdown("html/help_select.Rmd") ),
                        tabPanel("Help: Legend",  includeMarkdown("html/help_legend.Rmd") )
                        #
                      )
                    )

                    
                    
                  )
))


