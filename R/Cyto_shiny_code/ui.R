library(shiny)
#library(tidyr)
#this function takes r dataframes of nodes and edges and outputs the info/data in a list of 
# character strings which is then used by cytoscape.js to create the nework graph. 
source("./R_to_Cyjs_wrapper_functions/CCN2_pie.R")
# Cytoscape 'wrapper' function, uses output from CCN2_pie.R to create the cytoscape graph.
source("./R_to_Cyjs_wrapper_functions/cytoscapeJsSimpleNetwork3_strain.R")
# These are R functions used here in this shiny app but generic enough for non-shiny networking applications their uses 
# are found in "Create_Flu_nework#.R"    
source("./R_functions/FUNC_for_example_genomics3.R")
# all data is read in here.
source("model_data_read.R")
# This is the workhorse R script where the node and edge selections are made along with attributes, colors, shapes. Output from here is 
# used by CCN2_pie.R.
source("./R_functions/Create_Flu_network6.R")

# These will be needed when subnetworking is done.  
#library(BioNet)
#library(DLBCL)
#library(igraph)

#(theme = "bootstrap.css",
shinyUI(fluidPage(
        tags$head(
                #  Any javascript sourced in here must be placed in the "www" subdirectory in order for it to be found
                #  remember to add commas if more than one "tag", last one has no comma!
                tags$script(src="cytoscape.min.js")
                #  tags$script(src = 'http://cytoscape.github.io/cytoscape.js/api/cytoscape.js-latest/cola.v3.min.js'), 
        ),
        
        # Application title
        titlePanel("Fluomics RNA-seq Mouse Lung infected with H5N1/H1N1/H3N2"),
        
        sidebarLayout(
                sidebarPanel(width=3, 
                             sliderInput("q_value", "primary QVAL:", min=.001, max=.1, value=.001),
                             sliderInput("log2FC", "log 2 FC:", min=1.5 , max=3.0, value=2.0),
                             h4("Select virus(es)", align = "left"),
                             # The 'div(style="display:inline-block",)' portion simply allows the checkboxes to be displayed side by side, otherwise a new line will start for
                             # each line here.  
                             div( style = "color:#FF0000",checkboxInput("H1N1_s", "H1N1", value = T)),
                             
                             div(style="color:#0000FF", checkboxInput("H3N2_s", "H3N2", value = FALSE)),
                             
                             div(style="color:#006600", checkboxInput("H5N1_s", "H5N1", value = FALSE)),
                             
                             h5("Select time point(s)", align = "left"),
                             h6("Lighter color node = later time", align = "left"),
                             
                             div(style="display:inline-block", checkboxInput("T12H", "12H", value = FALSE)),
                             div(style="display:inline-block",checkboxInput("T24H", "24H", value = T)),
                             div(style="display:inline-block",checkboxInput("T48H", "48H", value = FALSE)),
                             div(style="display:inline-block",checkboxInput("T72H", "72H", value = FALSE)),
                             div(style="display:inline-block",checkboxInput("T96H", "96H", value = FALSE)),
                             
                             div(style="display:inline-block",checkboxInput("RN_ne", "Remove Nodes without Edges", value = TRUE)),
                             selectInput("layout", "Layout:", 
                                         #                         choices=c("grid", "random", "circle", "breadthfirst", "cose",
                                         #                                   "springy", "cola", "dagre", "force-directed"), 
                                         choices=c("grid", "random", "circle", "cose"), 
                                         
                                         selected="cose", 
                                         multiple=FALSE),
                             submitButton("Update Network"),
                             #            checkboxInput("addLinks", "Add Links on Nodes?", TRUE)
                             h5("Edge Color Legend", align = "left"),
                             
#                              node.colors.18 <- c("#CC3333", "#FF6600", "#FFFF33", "#00FF00", "#00CCCC","#00CCFF" ,"#3366FF",
#                                                  "#9933FF", "#FF00FF", "#FFCCCC", "#FFCC99", "#FFFFCC", "#CCFFCC", "#99FFCC",
#                                                  "#CCFFFF" ,"#99CCFF", "#CCCCFF", "#FFCCFF" )                             
                             # textOutput('color1'),
                             # textOutput('color2'),
                             #       This is the makeshift edge legend,  the outputs 'text1', 'text2' etc are the type of edge which is obtained from the currently displayed network 
                             #       (in server.R).   Next up is to hard assign the colors to edge type with check boxes.  
                             div( style = "color:#CC3333",checkboxInput("lit", "literature", value = T)),
                             div(style="color:#FF6600", checkboxInput("met", "metabolic", value = T)),
                             div(style="color:#FFFF33", checkboxInput("bin", "binary", value = T)),
                             div(style="color:#00FF00", checkboxInput("reg", "regulatory", value = T)),
                             div(style="color:#00CCCC", checkboxInput("com", "complexes", value = T)),
                             div(style="color:#00CCFF", checkboxInput("kin", "kinase", value = T)),
                             div(style="color:#3366FF", checkboxInput("sig", "signaling", value = T)),
                             #?checkboxInput
                             textOutput('text1'),
#                              tags$head(tags$style("#text1{color: #CC3333;
#                                           font-size: 20+px;
#                                           font-style: italic;
#                                           }"
#                              )
#                              ),
                             textOutput('text2'),
                             textOutput('text3'),
                             textOutput('text4'),
                             textOutput('text5'),
                             textOutput('text6'),
                             textOutput('text7')
                             
                ),
                
                # tags$head(tags$style("#text1{color: as.character(textOutput('color1'))  ;
                # generate ui elements from server.r, so that it can pick up names, colors and of course a varying list...
                mainPanel(div(id="cy"), htmlOutput("cytoscapeJsPlot"))
                
        )
))


