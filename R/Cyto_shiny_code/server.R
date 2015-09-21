library(shiny)


shinyServer(function(input, output, session) {
  # observe({ 
  output$ui1 <-  renderUI({
    if (is.null(input$selectSPECIES))
      return()
    
    # Depending on input$input_type, we'll generate a different
    # UI component and send it to the client.
    switch(input$selectSPECIES,
           "human" =   selectInput("selectDATA", label = h3("Data"), 
                                   choices = c("Human Genes HTBE" = 'human_rnaseq_HTBE', 
                                               "Human Genes MDM" = 'human_rnaseq_MDM'),
                                   #           "Human Proteomics ph" = 'proteomics_human_ph', 
                                   #            "Human Proteomics ub" = 'proteomics_human_ub'), 
                                   selected = 'human_rnaseq_HTBE'),
           "mouse" =   selectInput("selectDATA", label = h3("Data"), 
                                   choices = c("Mouse Genes" = 'mouse_rnaseq'), 
                                   #       "Mouse Proteomics ph" = 'proteomics_mouse_ph', 
                                   #        "Mouse Proteomics ub" = 'proteomics_mouse_ub', 
                                   #        "Mouse Metabolomics" = 'metabolomics_mouse'), 
                                   selected = 'mouse_rnaseq')
           
    )  # close switch
  }) # close output$ui
  #  }) # close observe
  
  #############
  ############
  DF1list <- reactive({ 
    if (is.null(input$selectDATA)) return('mouse_rnaseq')
    
    switch(input$selectDATA, 
           'human_rnaseq_HTBE'={
             withProgress(message = "loading human_rnaseq_data", value = 0.1, {
               load("./Linked_subdirectories/data/human_rnaseq_data")
             })
             data <- get_time(genomics_human)
             data <- data[data$cell_line == 'HTBE', ]
             MAP <- Human_map
             title <- "Fluomics RNA-seq Human HTBE Cell infected with H5N1/H1N1/H3N2"
             identity_attribute <- 'ellipse'
             DF1list <- namedList(data, MAP, title, identity_attribute)
           },
           'human_rnaseq_MDM'={
             
             withProgress(message = "loading human_rnaseq_data", value = 0.1, {
               load("./Linked_subdirectories/data/human_rnaseq_data")
             })
             data <- get_time(genomics_human)
             data <- data[data$cell_line == 'MDM', ]
             MAP <- Human_map
             title <- "Fluomics RNA-seq Human HTBE Cell infected with H5N1/H1N1/H3N2"
             identity_attribute <- 'ellipse'
             DF1list <- namedList(data, MAP, title, identity_attribute)
           },
           'mouse_rnaseq'={
             # case 'mouse_rnaseq' here...
             if(exists('human_rnaseq_data'))rm(human_rnaseq_data)
             #                     if(exists('Human_map'))rm(Human_map)
             withProgress(message = "loading mouse_rnaseq_data", value = 0.1, {
               load("./Linked_subdirectories/data/mouse_data_large3")
             })
             data <- get_time(genomics_mouse)
             MAP <- Mouse_map
             title <- "Fluomics RNA-seq Mouse Lung infected with H5N1/H1N1/H3N2"
             identity_attribute <- 'pentagon'
             DF1list <- namedList(data, MAP, title, identity_attribute)
           },
           DF1list # <- namedList(data, MAP, title)
    )
  })
  ############
  Time_Points <-  reactive({ sort(unique(DF1list()[["data"]]$time)) })
  #      Time_Points <-   sort(unique(DF1list()[["data"]]$time)) 
  output$TP1 <-  renderText({ Time_Points()[1] }) 
  output$TP2 <-  renderText({ Time_Points()[2] }) 
  output$TP3 <-  renderText({ Time_Points()[3] }) 
  output$TP4 <-  renderText({ Time_Points()[4] }) 
  output$TP5 <-  renderText({ Time_Points()[5] }) 
  output$title <-  renderText({ DF1list()[['title']] }) 
  
  DF2 <-  reactive({
    
    create.network1.1(DF1 = DF1list()[["data"]], Qval = input$q_value,  LOG2FC= input$log2FC)
    
    
  })
  
  network_object_react <- 
    reactive({  withProgress(message = "Selecting Nodes", value = 0.1, { create.network2.2 (DF3 = DF2()[['DF_reshape']], # Qval = input$q_value,  LOG2FC= input$log2FC,
                                                                                            h1n1s <- input$H1N1_s, h3n2s <- input$H3N2_s, h5n1s <- input$H5N1_s, map = DF1list()[["MAP"]], T1 = input$T1, 
                                                                                            T2 = input$T2, T3 = input$T3, T4 = input$T4, 
                                                                                            T5 = input$T5, time_points = Time_Points(), Show_other_viruses = input$SOV)
      # network_object_react <- create.network2(DF3 = DF2,  h1n1s <- T, h3n2s <- T, h5n1s <- T, t12h = T, t24h = T, t48h = F, t72h = F, t96h = F)
    }) # close 'with progress'
    })
  
  
  #        setProgress(value = 0.7,message = 'virus and time selection completed') 
  selected_edges <- reactive({ withProgress(message = "Selecting Edges", value = 0.6, {select_edges(edges = network_object_react()[[2]], lit <- input$lit, met <- input$met, bin <- input$bin, reg <- input$reg,
                                                                                                    com <- input$com, kin <- input$kin, sig <- input$sig)
    
    #                 selected_edges <- reactive({ select_edges(edges = edges, lit <- input$lit, met <- input$met, bin <- input$bin, reg <- input$reg,
    #                                          com <- input$com, kin <- input$kin, sig <- input$sig)
  }) # close 'with progress - selecting edges'
  }) # close reactive for edge selection
  
  
  nodeData <- reactive({ withProgress(message = "Removing Free Nodes", value = 0.7, {remove_free_nodes(rN_ne = input$RN_ne, nodes = network_object_react()[['AllNodes']] , edges = selected_edges(),
                                                                                                       node_id_colname = 'id')
  }) # close 'with progress - rfn'
  })
  ns2T <- reactive({ withProgress(message = "Removing Free Nodes", value = 0.7, {remove_free_nodes(rN_ne = input$RN_ne, nodes = network_object_react()[['ns2_display']] , edges = selected_edges(),
                                                                                                   node_id_colname = "id")
  }) # close 'with progress - rfn'
  })
  
  
  
  output$litL <- reactive({
    "literature" %in% network_object_react()[["edgeLegendData"]]$tags_for_colors               
  })
  output$metL <- reactive({
    "metabolic" %in% network_object_react()[["edgeLegendData"]]$tags_for_colors               
  })
  output$binL <- reactive({
    "binary" %in% network_object_react()[["edgeLegendData"]]$tags_for_colors               
  })
  output$comL <- reactive({
    "complexes" %in% network_object_react()[["edgeLegendData"]]$tags_for_colors               
  })
  output$kinL <- reactive({
    "kinase" %in% network_object_react()[["edgeLegendData"]]$tags_for_colors               
  })
  output$sigL <- reactive({
    "signaling" %in% network_object_react()[["edgeLegendData"]]$tags_for_colors               
  })
  output$regL <- reactive({
    "regulatory" %in% network_object_react()[["edgeLegendData"]]$tags_for_colors               
  })
  outputOptions(output, "litL", suspendWhenHidden=FALSE)
  outputOptions(output, "metL", suspendWhenHidden=FALSE)
  outputOptions(output, "binL", suspendWhenHidden=FALSE)
  outputOptions(output, "regL", suspendWhenHidden=FALSE)
  outputOptions(output, "comL", suspendWhenHidden=FALSE)
  outputOptions(output, "kinL", suspendWhenHidden=FALSE)
  outputOptions(output, "sigL", suspendWhenHidden=FALSE)
  outputOptions(output, "TP1", suspendWhenHidden=FALSE)
  
  #        output_x4 <- list(1,2)
  output$cytoscapeJsTable_edges <- DT::renderDataTable({ selected_edges() })
  output$cytoscapeJsTable_nodes <- DT::renderDataTable({ ns2T() })
  output$cytoscapeJsTable_nodes_attributes <- 
    DT::renderDataTable({ datatable(nodeData() )   })
  #                        output_x4 <- c(1,2)
  ### higlighted selections view
  #        observeEvent(input$cytoscapeJsTable_nodes_attributes_rows_selected, {
  #           output_x4 <- input$cytoscapeJsTable_nodes_attributes_rows_selected 
  #        })
  #        ns3T <- nodeData
  
  #         ns3T <- nodeData()[output_x4, ] 
  #         output$cytoscapeJsTable_nodes2 <- DT::renderDataTable({ ns3T  }) 
  #          nodeData_highlight <- nodeData() 
  #          nodeData_highlight[output_x4, "bkgrnd_highlight"] <- 'yellow'
  #          output$cytoscapeJsPlot_2 <-  renderPrint({
  #                 cyNetwork_h <- createCytoscapeNetwork(nodeData_highlight, selected_edges())
  #                 cytoscapeJsSimpleNetwork2(cyNetwork_h$nodes, cyNetwork_h$edges, layout=input$layout)
  #          }) # close 'render print'
  
  
  #                     }) # 'observe event'
  ##########
  withProgress(message = "Calculating", value = 0.1, {
    # nodeData <- nodeData_highlight
    output$cytoscapeJsPlot <-  renderPrint({ withProgress(message = "Generating Cytoscape.js", value = 0.8, {
      cyNetwork <- reactive({ createCytoscapeNetwork(nodeData(), selected_edges()) })
      cytoscapeJsSimpleNetwork2(cyNetwork()$nodes, cyNetwork()$edges, layout=input$layout)
    }) # close withProgress Generating Cytoscape.js
    }) 
  })  ## close 'with progress'
  
  
  
  
  
  
  
  #                 output$cytoscapeJsTable_nodes <- DT::renderDataTable({ ns2T() })
  #                 output$cytoscapeJsTable_edges <- DT::renderDataTable({ selected_edges() })
  
  output$cytoscapeJsTable_edges <- DT::renderDataTable({ selected_edges() })
  output$cytoscapeJsTable_nodes <- DT::renderDataTable({ ns2T() })
  output$cytoscapeJsTable_nodes_attributes <- 
    DT::renderDataTable({ datatable(nodeData() )   })
  
  
  #                  DT::renderDataTable({ datatable(nodeData(), options = list( selected = list(list(1,4))) )   })
  #                   DT::renderDataTable( { datatable( nodeData() ,
  #                                   selection = list(mode = 'multiple', selected = rownames(nodeData())[c(1,2)]) )  } )
  
  
  ###
  # print the selected indices
  #                 output$x4 = renderPrint({
  #                   s = input$cytoscapeJsTable_nodes_rows_selected
  #                   if (length(s)) {
  #                     cat('These rows were selected:\n\n')
  #                     cat(s, sep = ', ')
  #                   }
  #                 })
  
  
  
  #                   output_x4 <- reactive({ input$cytoscapeJsTable_nodes_rows_selected })
  
  #                  if(length(output_x4())){
  #                 ns3T <- reactive({ ns2T()[output_x4(), , drop = FALSE] })
  #                    output$cytoscapeJsTable_nodes2 <- DT::renderDataTable({ ns3T()  }) 
  #                  }
  
  
  #    output$cytoscapeJsTable_nodes2 <-render({ reactive({  ns2T()[output_x4 ,] }) }) 
  ###
  
  
  #          }
  
  
})


# head(cyNetwork[[1]], 1)
