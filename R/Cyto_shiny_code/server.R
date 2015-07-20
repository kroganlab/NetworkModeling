library(shiny)

shinyServer(function(input, output, session) {
        
        output$cytoscapeJsPlot <- renderPrint({
                
                DF2 <- reactive({
                        create.network1(DF1 = genomics_mouse, Qval = input$q_value,  LOG2FC= input$log2FC)
                })
                #  The following line is used for developing, debugging keeping everything in an R session without running the app.
                #  DATA2 <- create.network1(DATA = genomics_mouse, Qval = .001,  LOG2FC= 2.0)
                #  DF2 is DF2() simply because it is obtained from the 'reactive' expression above.
                network_object_react <- reactive({ create.network2(DF3 = DF2(), 
                                                                   h1n1s <- input$H1N1_s, h3n2s <- input$H3N2_s, h5n1s <- input$H5N1_s, t12h = input$T12H, 
                                                                   t24h = input$T24H, t48h = input$T48H, t72h = input$T72H, 
                                                                   t96h = input$T96H, rN_ne = input$RN_ne)
                                                   # reactive curly brace and parantheses for close below
                })
                
                
                #    The following line used for debugging / running in R.        
                #         network_object <- create.network2(DATA = DATA2,  h1n1s <- T,h3n2s <- F, h5n1s <- F,  t12h = T , t24h = T, t48h = F, t72h = F, t96h = F)        
                #  The following line is there only because network_object_react()   
                
                network_object <- network_object_react()
                nodeData <- network_object[["nodeData"]]
                edgeData <- network_object[[2]]
                
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
                if(num_edge_types > 7){output$text8 <- renderText({ as.character(network_object[[3]][8,2]) }) }
                if(num_edge_types > 8){output$text9 <- renderText({ as.character(network_object[[3]][9,2]) }) }
                if(num_edge_types > 9){output$text10 <- renderText({ as.character(network_object[[3]][10,2]) }) }
                if(num_edge_types > 10){output$text11 <- renderText({ as.character(network_object[[3]][11,2]) }) }
                
                # nodeData and edgeData are R dataframes with attributes (columns), createCytoscapeNetwork takes them and outputs 
                #   The nework as a list of character strings used by cytoscape.js (via the wrapper cytoscapeJsSimpleNetwork2)
                
                cyNetwork <- createCytoscapeNetwork(nodeData, edgeData)
                cytoscapeJsSimpleNetwork2(cyNetwork$nodes, cyNetwork$edges, layout=input$layout)
        }) 
        
})


# head(cyNetwork[[1]], 1)
