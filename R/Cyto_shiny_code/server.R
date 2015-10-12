library(shiny)


shinyServer(function(input, output, session) {
  # observe({ 
  output$ui1 <-  renderUI({
    if (is.null(input$selectSPECIES))
      return()
    # Depending on input$input_type, we'll generate a different
    # UI component and send it to the client.
    switch(input$selectSPECIES,
           "human" =   selectInput("selectCELLTYPE", label = h6("Cell Type"), 
                                   choices = c("HTBE" = 'HTBE', 
                                               "MDM" = 'MDM',
                                              "A549" = 'A549'), 
                                   selected = 'HTBE', multiple = F, selectize = F),
           "mouse" =   selectInput("selectCELLTYPE", label = h6("Cell Type"), 
                                   choices = c("Mouse Lung" = 'MOUSE_LUNG') , 
 #                                         "???" = 'mouse_lung1'), 
                                   selected = 'MOUSE_LUNG' , multiple = F, selectize = F)
    )  # close switch
  }) # close output$ui

  DF1list <- reactive({  
    dataname <- paste0(input$selectOMICS, "_", input$selectSPECIES)
    data_list <- list(); identity_attribute_list <- list()
    for(i in 1: length(dataname)){
      data1 <-  get_time(get(dataname[i]))
      data1 <- data1[data1$cell_line == input$selectCELLTYPE, ]
      data_list[[i]] <- data1
      identity_attribute_list[[i]] <- input$selectOMICS[i]
    }
  #  data <- data_list[[1]]
    MAP <- get(paste0(input$selectSPECIES, "_map"))
     DF1list <- namedList(data_list, MAP, identity_attribute_list)
  })
  
  ############
 # Time_Points <-  reactive({ sort(unique(DF1list()[["data"]]$time)) })
  Time_Points <-  reactive({ sort(unique(DF1list()[["data_list"]][[1]]$time)) })
  #      Time_Points <-   sort(unique(DF1list()[["data"]]$time)) 
  output$TP1 <-  renderText({ Time_Points()[1] }) 
  output$TP2 <-  renderText({ Time_Points()[2] }) 
  output$TP3 <-  renderText({ Time_Points()[3] }) 
  output$TP4 <-  renderText({ Time_Points()[4] }) 
  output$TP5 <-  renderText({ Time_Points()[5] }) 
  # output$title <-  renderText({ DF1list()[['title']] }) 

  
  DF2 <-  reactive({
    #  create.network1.1(DF1 = DF1list()[["data"]], Qval = input$q_value,  LOG2FC= input$log2FC)
    DF_list <- list()
    for(i in 1:length(DF1list()[["data_list"]])){
      DF_list[[i]] <- create.network1.1(DF1 = DF1list()[["data_list"]][[i]], Qval = input$q_value,  LOG2FC= input$log2FC)
    }
    DF_list
  })
  
  network_object_react <- 
    reactive({  withProgress(message = "Selecting Nodes", value = 0.1, { 
           N_O_list <- list()
            for(i in 1:length(DF2())){
                N_O_list[[i]] <-  create.network2.2 (DF3 = DF2()[[i]][['DF_reshape']], # Qval = input$q_value,  LOG2FC= input$log2FC,
                                 h1n1s <- input$H1N1_s, h3n2s <- input$H3N2_s, h5n1s <- input$H5N1_s, map = DF1list()[["MAP"]], T1 = input$T1, 
                                 T2 = input$T2, T3 = input$T3, T4 = input$T4, 
                                 T5 = input$T5, time_points = Time_Points(), Show_other_viruses = input$SOV, idmult = i, 
                                 shape_attribute = DF1list()[["identity_attribute_list"]][[i]] )
              
            }
           N_O_list
    }) # close 'with progress'
    })
  
  
  selected_edges <- reactive({ withProgress(message = "Selecting Edges", value = 0.6, {

             select_edges(edges = do.call("rbind",lapply(network_object_react(),function(x) x[["edgeData"]])), 
                          lit <- input$lit, met <- input$met, bin <- input$bin, reg <- input$reg,
                          com <- input$com, kin <- input$kin, sig <- input$sig)
             
                 
  }) # close 'with progress - selecting edges'
  }) # close reactive for edge selection
  
  
  nodeData <- reactive({ withProgress(message = "Removing Free Nodes", value = 0.7, 
                                    #  {remove_free_nodes(rN_ne = input$RN_ne, nodes = network_object_react()[[1]][['AllNodes']] ,
                                    {remove_free_nodes(rN_ne = input$RN_ne, nodes =  
                                                         do.call("rbind",lapply(network_object_react(),function(x) x[["AllNodes"]])),
                                                         edges = selected_edges(),  node_id_colname = 'id')
  }) # close 'with progress - rfn'
  })
  ns2T <- reactive({ withProgress(message = "Removing Free Nodes", value = 0.7, {remove_free_nodes(rN_ne = input$RN_ne, 
                                    nodes = do.call("rbind",lapply(network_object_react(),function(x) x[["ns2_display"]])), 
                                    edges = selected_edges(),  node_id_colname = "id")
  }) # close 'with progress - rfn'
  })
             
  # initial error might be here in the sequence of things, maybe edge selectors (widgets) need to be fully prepopulated in the Ui in order for the 
  # networking functions to grab the logic..  so maybe there is a 'null' loop then it keeps going and has what it needs?  Maybe use an observe statement somehow?
  # prob work on that "Coerced" error though since there is a chance that one could be having an effect.
  # observe({ print(network_object_react()[[1]][["edgeLegendData"]]) })
   # browser()
  
 # combined_edges <- reactive({ 
    Omics_edges <- reactive({ edges_between_omics(Nodes =  nodeData()) })
 #   combined_edges1 <-   reactive({ rbind(selected_edges(), Omics_edges()) })
 #   LENGTH_OMICS <- reactive({length(Omics_edges())})
#    if(LENGTH_OMICS == 0){
#       combined_edges <- reactive({ selected_edges() })
#    }
    combined_edges <- reactive({ assemble_edges(selected_edges(), Omics_edges()) })
    
 # })
  
#   if(length(Omics_edges() > 0)){
#     output <- combined_edges()
#   }else{
#     output <- selected_edges() 
#   }
#   output
 
 #  combined_edges <- reactive({ selected_edges()})   
      edge_types <- reactive({ do.call("rbind",lapply(network_object_react(),function(x) x[["edgeLegendData"]])) })
  output$litL <- reactive({
    "literature" %in% edge_types()$tags_for_colors 
  })
  output$metL <- reactive({
    "metabolic" %in% edge_types()$tags_for_colors            
  })
  output$binL <- reactive({
    "binary" %in% edge_types()$tags_for_colors               
  })
  output$comL <- reactive({
    "complexes" %in% edge_types()$tags_for_colors                
  })
  output$kinL <- reactive({
    "kinase" %in% edge_types()$tags_for_colors                
  })
  output$sigL <- reactive({
    "signaling" %in% edge_types()$tags_for_colors                
  })
  output$regL <- reactive({
    "regulatory" %in% edge_types()$tags_for_colors                
  })
  outputOptions(output, "litL", suspendWhenHidden=FALSE)
  outputOptions(output, "metL", suspendWhenHidden=FALSE)
  outputOptions(output, "binL", suspendWhenHidden=FALSE)
  outputOptions(output, "regL", suspendWhenHidden=FALSE)
  outputOptions(output, "comL", suspendWhenHidden=FALSE)
  outputOptions(output, "kinL", suspendWhenHidden=FALSE)
  outputOptions(output, "sigL", suspendWhenHidden=FALSE)
#  outputOptions(output, "TP1", suspendWhenHidden=FALSE)
  ## here we'll add the interomics edges just before rendering
  

  
  
  
  #########################
  #        output_x4 <- list(1,2)
#  output$cytoscapeJsTable_edges <- DT::renderDataTable({ selected_edges() })
#  output$cytoscapeJsTable_edges <- DT::renderDataTable({ combined_edges() })
#  output$cytoscapeJsTable_nodes <- DT::renderDataTable({ ns2T() })
#  output$cytoscapeJsTable_nodes_attributes <- 
#    DT::renderDataTable({ datatable(nodeData() )   })
 
  ##########
  withProgress(message = "Calculating", value = 0.1, {
  
    output$cytoscapeJsTable_nodes <-  DT::renderDataTable({ ns2T() }) 
    output$cytoscapeJsTable_nodes_attributes <- 
       DT::renderDataTable({ datatable(nodeData() )   }) 
 

    output$cytoscapeJsPlot <-   renderPrint({ withProgress(message = "Generating Cytoscape.js", value = 0.8, {

      cyNetwork <- createCytoscapeNetwork(nodeData(), combined_edges())
#      cyNetwork <- createCytoscapeNetwork(nodeData(), selected_edges())

      cytoscapeJsSimpleNetwork2(cyNetwork$nodes, cyNetwork$edges, layout=input$layout)
    }) # close withProgress Generating Cytoscape.js
    }) 
    output$cytoscapeJsTable_edges <-    ({ DT::renderDataTable({ combined_edges() }) }) 
#    }) # close cytoscapeJSplot 'reactive'
  })  ## close 'with progress'
  
  
 


  
  
})


# head(cyNetwork[[1]], 1)
