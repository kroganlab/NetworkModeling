#target_time <- "01D"; strain = "all"
# create.network (used to be create.network1) is now only used for reshape BEFORE shiny app, so only once on a new data set.
get_time <- function(DF1 = genomics_mouse){
  #        load("./Linked_subdirectories/data/human_rnaseq_data"); DF1 <- get_time(genomics_human);      map <- Human_map
  #       load("./Linked_subdirectories/data/mouse_data_large3"); DF1 <- get_time(genomics_mouse);      map <- Mouse_map
  gm_subf <- DF1 # data.table(DF1)
  # set up time column, criteria, and remove duplicate rows.
  gm_subf$time <-  gsub("(.*)_(.*)" , '\\2', gm_subf$condition_2 )
  gm_subf[gm_subf$time == '01D', "time"] <- "24H"; gm_subf[gm_subf$time == '02D', "time"] <- "48H"
  gm_subf[gm_subf$time == '03D', "time"] <- "72H"; gm_subf[gm_subf$time == '04D', "time"] <- "96H"
  # time_points <- sort(unique(gm_subf$time))
  gm_subf
}


create.network1.1 <- function( DF1 = genomics_mouse,   Qval = .001,  LOG2FC = 2.0)
{
  # ; Qval = .1;  LOG2FC = 1.5
  gm_subf <- DF1  # not if running manually.
  gm_subf <- unique(gm_subf)
  gm_subf <- gm_subf[!(is.na(gm_subf$entrez_id) == T & is.na(gm_subf$uniprot_ac) == T),]
  gm_subf <- gm_subf[((gm_subf$q_value <= Qval) & ((gm_subf$log2fc <= -1*LOG2FC) | (gm_subf$log2fc >= LOG2FC ))), ]
  #       CELL_LINE <- 'MDM' # 'A549, HTBE, MDM'
  #       gm_subf <- gm_subf[gm_subf$cell_line == CELL_LINE, ]
  gm_subf <- unique(gm_subf)
  
  # the following function is used by the aggregate command to find the worst fold change for entries from 2 or more experiments.
  
        min_abs2 <- function(Y){
          val <-  Y[which.max( abs(Y) )]
          return(val)
        }
  
  
  
  M <- aggregate(log2fc ~ entrez_id + strain + time + symbol, data = gm_subf[ , c('entrez_id', 'strain', 'time',  'symbol', 'log2fc')], min_abs2 )
  
  # time_points <<- sort(unique(gm_subf$time))
  
  # rm(gm_subf)
  M <- unique(M); M <- M[M$strain != 'WSN' , ];  M <- data.table(M)
  M <- unique(M)
  

  
  
  M[is.na(M$log2fc) , 'log2fc' ] <- 0
  
  # M2E is created to make sure all 3 viruses are present in the data set before it is reshaped (dcast).  This way the next function will work properly.  These added rows don't
  # show up in the final plot or data tables - they just simply keep everything working properly.
  M2E <- M[1:3, ]; 
  M2E[1, 'entrez_id'] <- 0; M2E[1, 'symbol'] <- 'fill - ignore';  M2E[1, 'strain'] <- 'H1N1'; M2E[1, 'log2fc'] <- 0.0
  M2E[2, 'entrez_id'] <- 0; M2E[2, 'symbol'] <- 'fill - ignore';  M2E[2, 'strain'] <- 'H3N2';  M2E[2, 'log2fc'] <- 0.0
  M2E[3, 'entrez_id'] <- 0; M2E[3, 'symbol'] <- 'fill - ignore';  M2E[3, 'strain'] <- 'H5N1';  M2E[3, 'log2fc'] <- 0.0
  M <- rbind(M, M2E)
  
  DF_reshape <- dcast(M, entrez_id + symbol + time ~ strain, value.var = 'log2fc'  )
  namedList(DF_reshape)
}



# this is the main workhorse function.  time and virus selections are made but most of the work is done for the pie color chart on the nodes.
create.network2.2 <- function( DF3 = DF_reshape,  #  Qval = .001,  LOG2FC = 2.0,
                               h1n1s = F, h3n2s = F, h5n1s = F, map = Mouse_map,
                               T1 = F, T2 = F, T3 = F, T4 = F, T5 = F, T6 = F, T7 = F, T8 = F, time_points,
                               Show_other_viruses = F , idmult = 1, shape_attribute = 'genomics') #Node_stats = NULL)
{
  # source in the following line for running directly.
  # DF3 = DF_reshape; h1n1s = T; h3n2s = T; h5n1s = T; T1 = F; T2 = T; T3 = T; T4 = F; T5 = T;  Show_other_viruses = T
  
  t1 <- 9; t2 <- 9; t3 <- 9; t4 <- 9; t5 <- 9; t6 <- 9
  h1n1 <- 9; h3n2 <- 9 ; h5n1 <- 9
  if(h1n1s == T){ h1n1 <- 1}; if(h3n2s == T){ h3n2 <- 1} ; if(h5n1s == T){ h5n1 <- 1} 
  #     time_points <<- sort(unique(DF1$time))
  if(T1 == T){ t1 <- time_points[1]}; if(T2 == T){ t2 <- time_points[2]}
  if(T3 == T){ t3 <- time_points[3]}; if(T4 == T){ t4 <- time_points[4]}
  if(T5 == T){ t5 <- time_points[5]}
  DF <- DF3  
  max_time <- max(as.numeric(unique(substr(DF$time, 1,2))))
  
  
  # This set of if statements simply sets columns equal to zero for any virues not chosen if 'show other viruses' = F.
  if(Show_other_viruses == F){
    if(h1n1s == F){
      DF$H1N1 <- 0
    }
    if(h3n2s == F){
      DF$H3N2 <- 0
    }
    if(h5n1s == F){
      DF$H5N1 <- 0
    }
  }
  
  DF[is.na(DF$H1N1) , 'H1N1' ] <- 0; DF[is.na(DF$H5N1) , 'H5N1' ] <- 0; DF[is.na(DF$H3N2) , 'H3N2' ] <- 0
  
  # set up node data for cytoscape attributes (2 copies, one for display as tab, other used in cytoscape logic for colors etc)   
  
  DF_display <- DF # DF <- DF_display
  
  DF[DF$H1N1!=0 , 'H1N1' ] <- 1;  DF[DF$H5N1!=0 , 'H5N1' ] <- 1;  DF[ DF$H3N2!=0 , 'H3N2' ] <- 1
  selected_rows_by_strain <- DF$H1N1 == h1n1 | DF$H5N1 == h5n1 | DF$H3N2 == h3n2
  ns1 <- DF[selected_rows_by_strain, ]; ns1_display <- DF_display[selected_rows_by_strain, ] 
  selected_rows_by_time <- ns1$time == t1 | ns1$time == t2 | ns1$time == t3 | ns1$time == t4 | ns1$time == t5
  ns2 <- ns1[selected_rows_by_time, ] ;    ns2_display <- ns1_display[selected_rows_by_time, ]
  #  ns2_display will be for the display, ns2_displayL will be used ultimately to set up logic 
  # for attributes (size of node based on fold change) 
  
  ## logic to remove or keep 'unselected' virus strains that occur for the same genes at the same chosen time points. 
  ## The node is still there, (i.e no rows are removed)  
  
  ns3 <- ns2
  # rm(ns1,ns2) # save ns2 for output
  rm(ns1)
  names(ns2_display)[names(ns2_display)=='entrez_id'] <- "id"; names(ns2_display)[names(ns2_display)=='symbol'] <- "name"
  names(ns3)[names(ns3)=='entrez_id'] <- "id"; names(ns3)[names(ns3)=='symbol'] <- "name"
  
  # the following lines are setting up the columns 'h1n1_time' which will ultimately be used to shade the color of the node (or pie portion of the node)
  # for the appropriate time.  
 
  ns3$time <- as.numeric(substr(ns3$time, 1,2)) 
  ns3$h1n1_time <- (ns3$H1N1)*ns3$time; ns3$h3n2_time <- (ns3$H3N2)*ns3$time
  ns3$h5n1_time <- (ns3$H5N1)*ns3$time
  ns3$h1n1_time <- ifelse(ns3$h1n1_time == 0, 99, ns3$h1n1_time*100/max_time ) 
  ns3$h3n2_time <- ifelse(ns3$h3n2_time == 0, 99, ns3$h3n2_time*100/max_time ) 
  ns3$h5n1_time <- ifelse(ns3$h5n1_time == 0, 99, ns3$h5n1_time*100/max_time )
  
  ns3$name_id <- paste0(ns3$name, "_", ns3$id)
  ns3 <- ns3[,colnames(ns3[ , !names(ns3) %in% c("egg", "WSN", 'time')])] 
  # now for the display data:
  ns2_display$name_id <- paste0(ns2_display$name, "_", ns2_display$id)
  ns2_display3 <- ns2_display[,colnames(ns2_display[ , !names(ns2_display) %in% c("egg", "WSN", 'time')])] 
  
  # This line finds the earliest time point so that the shading of the color on the node reflects that earlier time point.
  ns3_h1n1min <- aggregate(h1n1_time ~ name_id, data = ns3, min )
  # the next two lines simply result in a '1' if something exists for that omics entity (name_id) for that viurs
  ns3_h1n1sum <- aggregate(H1N1 ~ name_id, data = ns3, sum )
  ns3_h1n1sum[ns3_h1n1sum$H1N1 > 0, "H1N1" ] = 1
  # the following line fills in the appropriate value for node color shading.
  ns3_h1n1_m_s <- merge(ns3_h1n1min, ns3_h1n1sum, by = 'name_id')
  
  # same thing for the remaining viruses.
  ns3_h3n2min <- aggregate(h3n2_time~ name_id, data = ns3, min )
  ns3_h3n2sum <- aggregate(H3N2 ~ name_id, data = ns3, sum )
  ns3_h3n2sum[ns3_h3n2sum$H3N2 > 0, "H3N2" ] = 1
  ns3_h3n2_m_s <- merge(ns3_h3n2min, ns3_h3n2sum, by = 'name_id')
  
  ns3_h5n1min <- aggregate(h5n1_time~ name_id, data = ns3, min )
  ns3_h5n1sum <- aggregate(H5N1 ~ name_id, data = ns3, sum )
  ns3_h5n1sum[ns3_h5n1sum$H5N1 > 0, "H5N1" ] = 1
  ns3_h5n1_m_s <- merge(ns3_h5n1min, ns3_h5n1sum, by = 'name_id')
  
  # set up some fold change attribute using stuff from display df.
  min_abs2 <- function(Y){
    val <-  Y[which.min( abs(Y) )]
    return(val)
  }
  min_abs2a <- function(Y){
    val <-  Y[which.min( abs(Y) )]
    return(abs(val))
  }
  
  min_abs2s <- function(Y){
    val <-  Y[which.max( abs(Y) )]
    return(sign(val))
  }
  
  ns2_displayL <- ns2_display3
  ns2_displayL[abs(ns2_displayL$H1N1) <= 0.01 , 'H1N1' ] <- 999;  ns2_displayL[abs(ns2_displayL$H5N1) <= 0.01 , 'H5N1' ] <- 999
  ns2_displayL[ abs(ns2_displayL$H3N2) <= 0.01 , 'H3N2' ] <- 999
  # this line takes the lowest fold change (abs value) from the three viruses
  ns2_displayL$log2fc_min_w_sign <- apply(ns2_displayL[ , c('H1N1', 'H3N2', 'H5N1')] ,  1, min_abs2 )
  # ns2_displayL$log2fc_sign <- apply(ns2_displayL[ , c('H1N1', 'H3N2', 'H5N1')] ,  1, min_abs2s )
  # Now its aggregated where there are multple time points (again min taken)
  # so instead do 2 aggregates, one for the min(do abs now) and one for the sign, both will have to be merged in.
  ns2_displayL_aggregate <- aggregate(log2fc_min_w_sign ~ name_id, data = ns2_displayL, min_abs2a )
  ns2_displayL_aggregate_sign <- aggregate(log2fc_min_w_sign ~ name_id, data = ns2_displayL, min_abs2s )
  # now change colname back to meaningful one"
  names(ns2_displayL_aggregate)[names(ns2_displayL_aggregate) == 'log2fc_min_w_sign'] <- 'log2fc_attrib'
  # plop in colors for sign here!!
  ns2_displayL_aggregate_sign$sign_color <-  'white'
  ns2_displayL_aggregate_sign[ns2_displayL_aggregate_sign$log2fc_min_w_sign <= 0.0, "sign_color"] <- "orange"
  ns2_displayL_aggregate_sign <- ns2_displayL_aggregate_sign[ , c('name_id', 'sign_color' )]
  
  ns3a <- ns3[ , !names(ns3) %in% c("h1n1_time","h3n2_time","h5n1_time", 'H5N1', "H3N2", "H1N1", "WSN", 'time')]
  ns3a <- merge(ns3_h1n1_m_s, ns3a, by = 'name_id')
  ns3a <- merge(ns3_h3n2_m_s, ns3a, by = 'name_id')
  ns3a <- merge(ns3_h5n1_m_s, ns3a, by = 'name_id')
  ns3a <- unique(ns3a); ns2_display <- unique(ns2_display)
  # merge in (probably in same row order, but just in case) log2fc attribute:
  ns3b <- merge(ns3a, ns2_displayL_aggregate, by = 'name_id')
  ns3c <- merge(ns3b, ns2_displayL_aggregate_sign, by = 'name_id')
  ns4 <- ns3c; rm(ns3a, ns3b)
  
  
  ns4$sumvirus <-  rowSums(ns4[ , c('H1N1' , 'H3N2' , 'H5N1')]); ns5 <- ns4
  ns5$H1N1 <- (18/ns4$sumvirus)*ns4$H1N1; ns5$H5N1 <- (18/ns4$sumvirus)*ns4$H5N1; ns5$H3N2 <- (18/ns4$sumvirus)*ns4$H3N2
  #  ns5 <- unique(ns5[ , -2]); rm(ns4)
  ns5 <-  ns5[ , !names(ns5) %in% c("time","u")]; rm(ns4)
  
  # now put time as a decimal percent and reversed - 12 = 1.0 and 96 <- .25?
  ns6 <- ns5; 
  ns6$h1n1_time <- 1- ((ns5$h1n1_time - 12)/ 100); ns6$h3n2_time <- 1- ((ns5$h3n2_time - 12)/ 100)
  ns6$h5n1_time <- 1- ((ns5$h5n1_time - 12)/ 100) 
  #      ns6 <- ns6[, -9]; rm(ns5)
  ns6 <-  ns6[ , !names(ns6) %in% c("sumvirus","name_id")]; rm(ns5)
  AllNodes <- ns6
  
  AllNodes$bkgrnd_highlight <- 'white'
  
  edges <- Prey_connections_OR3(map, "geneid_1", "geneid_2", AllNodes , "id", sourcecol = "source") 
  edges <- edges[edges$geneid_1 != edges$geneid_2, ] # remove self edges.
  #  Get only those connections which are in the same time(s), same strain(s)...
  edges <-      edges[edges$geneid_2 %in% ns6$id, ]
  
  colnames(edges) <- c("source", "target", "info_source")
  edges <- edges[ , c("source", "target", "info_source") ]
  
  edges2 <- edges[ ,c('source', 'target')]
  edges3 <- apply(edges2, 1 , sort)
  
  edges3.T <- t(edges3[,1:ncol(edges3)])
  edges3 <- as.data.frame(edges3.T)
  colnames(edges3) <- c("source", 'target')
  edges3$info_source <- edges$info_source
  edgeData <- unique(edges3); rm(edges3)
  
  ### move Expand edges and remove unconected nodes outside of this function into new one, this way the logic on which edges can be toggled after the network
  ## has been initially selected.  (we could remove edge types from consideration upstream I suppose, but that might be more complicated in the end The idea 
  # here is that it is relatively fast when reselecting based on edge type )
  Expanded_edges <- Add_edges_from_mult_sources(edgeData, "info_source", split_string = ";") 
  rm(edgeData)
  #      edgeList<- color_column2(Expanded_edges, "info_source")
  edgeList<- color_column3(Expanded_edges, "info_source")
  
  edgeData <- edgeList[[1]]; edgeLegendData <- edgeList[[2]]
  edgeData <- unique(edgeData)
  multiplier <- as.numeric(idmult)
   edgeData$source <- edgeData$source + (1000000000*multiplier )
   edgeData$target <- edgeData$target + (1000000000*multiplier )
   edgeData$thick <- 1.0
   AllNodes$shape <- switch(shape_attribute, "genomics" = "ellipse", "proteomics" = "hexagon", "metabolomics" = "triangle")
   AllNodes$id <- as.numeric(AllNodes$id) + (1000000000*multiplier )
   ns2_display$id <- as.numeric(ns2_display$id) + (1000000000*multiplier) 
  edge_list <- namedList(AllNodes, edgeData, edgeLegendData, ns2_display )
  
  edge_list
}



edges_between_omics <- function(Nodes){
  Nodes$colA <- Nodes$id
  Nodes$colB <- Nodes$id - (floor(Nodes$id/1000000000) * 1000000000)
  library(plyr)
  #  colA <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n")
  #  colB <- c(1,2,3,1,1,4,4,8,3,1,3,5,6,7)
  #  DF <- data.frame(colA, colB)
  DF <- Nodes[ , c('colA', 'colB')]
  list_fa <- unique(DF[duplicated(DF$colB), 'colB' ])
  data1 <- list()
  if(length(list_fa) > 0){
    for(i in 1:length(list_fa)){
      combination <- combn(DF[DF$colB == list_fa[[i]], "colA"], m = 2)
      data1[[i]] <- data.frame(source = combination[1,], target = combination[2,])
    }
  }
  inter_omics_edges <- rbindlist(data1)
  if(length(list_fa) > 0){    
    inter_omics_edges <- as.data.frame((inter_omics_edges))
    inter_omics_edges$info_source <- "inter_omics"
    #  inter_omics_edges$color <- "#FF00FF"
    inter_omics_edges$color <- "#000000"
    inter_omics_edges$thick <- 5.0
    inter_omics_edges <- inter_omics_edges[inter_omics_edges$source != inter_omics_edges$target, ] # remove self edges.
    inter_omics_edges <- unique(inter_omics_edges)
  }
  inter_omics_edges
}

assemble_edges <- function(selected_edges, Omics_edges){
  if(length(Omics_edges) > 0){
    output <- rbind(selected_edges, Omics_edges)
  }else{
    output <- selected_edges
  }
  output
}

remove_free_nodes <- function(rN_ne, nodes, edges, node_id_colname){
  if(rN_ne){
    nodesout <- nodes[nodes[ , node_id_colname ] %in% as.vector(unlist(edges[ , c('source', 'target')])) , ]    
  }else{
    nodesout<- nodes
  }   
  nodesout
}
# edges <- edgeData
select_edges <- function(edges, lit = F, met = F, bin = F, reg = F, com = F, kin = F, sig = F){
  #        lit <- T; met <- T; bin <- T; reg <- T; com <- T; kin <- T; sig <- T        
  if(lit){lit <- 'literature'}; if(met){met <- 'metabolic'}; if(bin){bin <- "binary"}; if(reg){reg <- 'regulatory'}
  if(com){com <-'complexes'}; if(kin){kin <- 'kinase'}; if(sig){sig <- 'signaling'}
  string <-  paste0(lit,"|" ,met,"|" ,bin,"|" ,reg,"|" ,com, "|" ,kin,"|" ,sig)
  logic <- grepl(string, edges$info_source)
  edges <- edges[logic,]
  edge_key <- data.frame(info_source = c("literature","metabolic","binary","regulatory","complexes","kinase","signaling"),
                         color = c("#FF0000", "#000000", "#0000FF", "#008000", "#800000","#8D38C9" ,"#F88017"))
  edges <- merge(edges[ , c('source', 'target', 'info_source', 'thick')], edge_key, by.x = 'info_source')
  edges
}


namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}