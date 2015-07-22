library(shiny)

shinyServer(function(input, output, session) {
  #      ?withProgress
        output$cytoscapeJsPlot <- renderPrint({
                withProgress(message = "Selecting log2FC, Qval, performing reshape", value = 0.5, {
                        DF2 <- reactive({
                                create.network1(DF1 = genomics_mouse, Qval = input$q_value,  LOG2FC= input$log2FC)
                        })
                        DF2reactive <- DF2()
                }) # close "withProgress" : "Selecting log2FC, Qval, performing reshape".
                        #  The following line is used for developing, debugging keeping everything in an R session without running the app.
                        #  DF2 <- create.network1(DF1 = genomics_mouse, Qval = .001,  LOG2FC= 2.0)
                        #  DF2 is DF2() simply because it is obtained from the 'reactive' expression above.
                  
                withProgress(message = "Selecting Virus and timepoints", value = 0.5, {         
               network_object_react <- reactive({ create.network2(DF3 = DF2reactive, # Qval = input$q_value,  LOG2FC= input$log2FC,
                                                                   h1n1s <- input$H1N1_s, h3n2s <- input$H3N2_s, h5n1s <- input$H5N1_s, t12h = input$T12H, 
                                                                   t24h = input$T24H, t48h = input$T48H, t72h = input$T72H, 
                                                                   t96h = input$T96H)
                                                   # reactive curly brace and parantheses for close below
                })
                
                
                #    The following line used for debugging / running in R.        
                #         network_object <- create.network2(DF3 = DF2,  h1n1s <- T,h3n2s <- F, h5n1s <- F,  t12h = T , t24h = T, t48h = F, t72h = F, t96h = F)        
                #  The following line is there only because network_object_react()   
                
                network_object <- network_object_react()
                # network_object <- network_list
               
                }) # close "Selecting Virus and timepoints"
                allnodes <- network_object[['AllNodes']]
                edges <- network_object[[2]]
                
                
                selected_edges <- reactive({ select_edges(edges = edges, lit <- input$lit, met <- input$met, bin <- input$bin, reg <- input$reg,
                                         com <- input$com, kin <- input$kin, sig <- input$sig)
                
                }) # close reactive for edge selection
                edgeData <- selected_edges()
                #edgeData <- select_edges(edges = edges, lit <- T, met <- T, bin <- T, reg <- T, com <- T, kin <- T, sig <- T)
                #edgeData <- select_edges(edges = edges, lit <- F, met <- F, bin <- F, reg <- F, com <- F, kin <- F, sig <- F)
                nodeData <- remove_free_nodes(rN_ne = input$RN_ne, nodes = allnodes , edges = edgeData)
                 # nodeData <- remove_free_nodes(rN_ne = T, nodes = allnodes , edges = edgeData)
                #  this is really messy here, basically just grabbing info for the makeshift edge legend, which should be updated to checkboxes soon 
                # meaning legend no longer dynamic, or possibly.. but still with checkboxes and assigned colors AND the option for edge thickness tied to 
                # # of associations.
                
                num_edge_types <- nrow(network_object[[3]])
                if(num_edge_types > 0){output$text1 <- renderText({ as.character(network_object[[3]][1,2]) }) }
                if(num_edge_types > 1){output$text2 <- renderText({ as.character(network_object[[3]][2,2]) }) }
                if(num_edge_types > 2){output$text3 <- renderText({ as.character(network_object[[3]][3,2]) }) }
                if(num_edge_types > 3){output$text4 <- renderText({ as.character(network_object[[3]][4,2]) }) }
                if(num_edge_types > 4){output$text5 <- renderText({ as.character(network_object[[3]][5,2]) }) }
                if(num_edge_types > 5){output$text6 <- renderText({ as.character(network_object[[3]][6,2]) }) }
                if(num_edge_types > 6){output$text7 <- renderText({ as.character(network_object[[3]][7,2]) }) }
               
                # nodeData and edgeData are R dataframes with attributes (columns), createCytoscapeNetwork takes them and outputs 
                #   The nework as a list of character strings used by cytoscape.js (via the wrapper cytoscapeJsSimpleNetwork2)
               withProgress(message = 'creating input for Cytoscape.JS', value = 0.5, {
                cyNetwork <- createCytoscapeNetwork(nodeData, edgeData)
               }) # close 'creating input for Cytoscape.JS'
               withProgress(message = 'creating Cytoscape.JS plot', value = 0.5, {
                cytoscapeJsSimpleNetwork2(cyNetwork$nodes, cyNetwork$edges, layout=input$layout)
              #  )
                         }) # for 'creating Cytoscape.JS plot'
  }     )    
})


# head(cyNetwork[[1]], 1)
