#target_time <- "01D"; strain = "all"
# create.network (used to be create.network1) is now only used for reshape BEFORE shiny app, so only once on a new data set.
create.network1 <- function( DF1 = genomics_mouse,   Qval = .001,  LOG2FC = 2.0)
{
        # DF1 = genomics_mouse; Qval = .001;  LOG2FC = 2.0
        gm_subf <- DF1
        # set up time column, criteria, and remove duplicate rows.
        gm_subf$time <- substr(gm_subf$condition_2,6,8)
        gm_subf[gm_subf$time == '01D', "time"] <- "24H"; gm_subf[gm_subf$time == '02D', "time"] <- "48H"
        gm_subf[gm_subf$time == '03D', "time"] <- "72H"; gm_subf[gm_subf$time == '04D', "time"] <- "96H"
        gm_subf$criteria <- 1
        gm_subf <- unique(gm_subf)
        gm_subf <- gm_subf[!(is.na(gm_subf$entrez_id) == T & is.na(gm_subf$uniprot_ac) == T),]
        gm_subf <- gm_subf[((gm_subf$q_value <= Qval) & ((gm_subf$log2fc <= -1*LOG2FC) | (gm_subf$log2fc >= LOG2FC ))), ]
        
        M <- gm_subf[ , c('entrez_id', 'strain', 'time', 'criteria', 'symbol' )] # c("strain", idtype, "criteria")
        rm(gm_subf)
        M <- unique(M)
        DF_reshape <- reshape(M, timevar = "strain" ,  idvar = c('entrez_id','symbol' , 'time'), direction = "wide")
        #         DF_reshape$criteria <- 1
        DF_reshape[is.na(DF_reshape)] <- 0
        #         # names = gsub("criteria.", "", colnames(DF_reshape[,c(2:ncol(DF_reshape))]))
        #         DF_reshape2 <- reshape(DF_reshape, timevar = "time" , 
        #                                idvar = c('entrez_id','symbol' ,'criteria.H5N1', 'criteria.H1N1', 'criteria.H3N2'), direction = "wide")
        #         DF_reshape2[is.na(DF_reshape2)] <- 0
        #         names = gsub("criteria.", "n", colnames(DF_reshape2[,c(1:ncol(DF_reshape2))]))
        #        colnames(DF_reshape2) <- names
        names = gsub("criteria.", "", colnames(DF_reshape[,c(1:ncol(DF_reshape))]))
        colnames(DF_reshape) <- names
        DF_reshape
}

 
create.network2 <- function( DF3 = DF_reshape,  #  Qval = .001,  LOG2FC = 2.0,
                             h1n1s = F, h3n2s = F, h5n1s = F,
                             t12h = F, t24h = F, t48h = F, t72h = F, t96h = F)
{
        
        # DF3 <- DF2 ;DF3 = DF_reshape; h1n1s = T; h3n2s = F; h5n1s = F; t12h = T; t24h = T; t48h = T; t72h = T; t96h = T ; Qval = .1001; LOG2FC = 2.0
        
        
        t12 <- 9; t24 <- 9; t48 <- 9; t72 <- 9; t96 <- 9
        h1n1 <- 9; h3n2 <- 9 ; h5n1 <- 9
        if(h1n1s == T){ h1n1 <- 1}; if(h3n2s == T){ h3n2 <- 1} ; if(h5n1s == T){ h5n1 <- 1} 
        if(t12h == T){ t12 <- '12H'}; if(t24h == T){ t24 <- '24H'}  ;if(t48h == T){ t48 <- '48H'}   
        if(t72h == T){ t72 <- '72H'}  ;if(t96h == T){ t96 <- '96H'}         
         DF <- DF3  
  #      DF <- DF3[((DF3$q_value <= Qval) & ((DF3$log2fc <= -1*LOG2FC) | (DF3$log2fc >= LOG2FC ))), ]
    # DF_view <- DF_reshape[2000:3000,]    
        ns1 <- DF[DF$H1N1 == h1n1 | DF$H5N1 == h5n1 | DF$H3N2 == h3n2, ]
        ns2 <- ns1[ns1$time == t12 | ns1$time == t24 | ns1$time == t48 | ns1$time == t72 | ns1$time == t96, ]
        ns3 <- ns2; rm(ns1,ns2)
        colnames(ns3)[1] <- "id"; colnames(ns3)[2] <- "time"; colnames(ns3)[3] <- "name"
        ns3$time <- as.numeric(substr(ns3$time, 1,2)) 
        ns3$h1n1_time <- (ns3$H1N1)*ns3$time; ns3$h3n2_time <- (ns3$H3N2)*ns3$time
        ns3$h5n1_time <- (ns3$H5N1)*ns3$time
        ns3$h1n1_time <- ifelse(ns3$h1n1_time == 0, 99, ns3$h1n1_time ) 
        ns3$h3n2_time <- ifelse(ns3$h3n2_time == 0, 99, ns3$h3n2_time ) 
        ns3$h5n1_time <- ifelse(ns3$h5n1_time == 0, 99, ns3$h5n1_time )
        #   ns3[ , c("id", "h1n1_time", "h3n2_time", "h5n1_time")]
        ns4 <- ns3; rm(ns3)
        #ns4[ns4$id == ns4[3,"id"],   ]; rows <- ns4$id == ns4[3,"id"]; min(ns4[rows, "h1n1_time"])
        for(i in 1:nrow(ns4)){
                if(nrow(ns4[ns4$id == ns4[i,"id"], ]) > 1){
                        rows <- ns4$id == ns4[i,"id"]
                        mh1n1t <- min(ns4[rows, "h1n1_time"])
                        ns4[rows, "h1n1_time"] <- mh1n1t
                        ns4[rows, "H1N1"] <- ifelse(sum(ns4[rows ,"H1N1" ]) > 0 , 1, 0)
                        mh3n2t <- min(ns4[rows, "h3n2_time"])
                        ns4[rows, "h1n1_time"] <- mh1n1t
                        ns4[rows, "H3N2"] <- ifelse(sum(ns4[rows ,"H3N2" ]) > 0 , 1, 0)
                        mh5n1t <- min(ns4[rows, "h5n1_time"])
                        ns4[rows, "h5n1_time"] <- mh5n1t
                        ns4[rows, "H5N1"] <- ifelse(sum(ns4[rows ,"H5N1" ]) > 0 , 1, 0)
                }
        }
        
        ns4$sumvirus <-  rowSums(ns4[ , c('H1N1' , 'H3N2' , 'H5N1')]); ns5 <- ns4
        ns5$H1N1 <- (18/ns4$sumvirus)*ns4$H1N1; ns5$H5N1 <- (18/ns4$sumvirus)*ns4$H5N1; ns5$H3N2 <- (18/ns4$sumvirus)*ns4$H3N2
        ns5 <- unique(ns5[ , -2]); rm(ns4)
        
        # now put time as a decimal percent and reversed - 12 = 1.0 and 96 <- .25?
        ns6 <- ns5; 
        ns6$h1n1_time <- 1- ((ns5$h1n1_time - 12)/ 100); ns6$h3n2_time <- 1- ((ns5$h3n2_time - 12)/ 100)
        ns6$h5n1_time <- 1- ((ns5$h5n1_time - 12)/ 100) 
        ns6 <- ns6[, -9]; rm(ns5)
        AllNodes <- ns6
        edges <- Prey_connections_OR3(Mouse_map, "geneid_1", "geneid_2", AllNodes , "id", sourcecol = "source") 
        edges <- edges[edges$geneid_1 != edges$geneid_2, ] # remove self edges.
        #  Get only those connections which are in the same time(s), same strain(s)...
        edges <-      edges[edges$geneid_2 %in% ns6$id, ]
        
        
        #        nodeData <- ns3[ ,c("id", "name", "time", "H1N1", "H3N2", "H5N1")]
        #        nodeData <- ns3
        colnames(edges) <- c("source", "target", "info_source")
        edges <- edges[ , c("source", "target", "info_source") ]

  #  Remove duplicate edges where node order is reversed:
        if(nrow(edges) > 0){
                edgeData <- remove_dup_edge_rev_order(edges_sort = edges, colname1 = 'source', colname2 = 'target', attribute_colname1 = "info_source")
                
        }
        ### move Expand edges and remove unconected nodes outside of this function into new one, this way the logic on which edges can be toggled after the network
        ## has been initially selected.  (we could remove edge types from consideration upstream I suppose, but that might be more complicated in the end The idea 
        # here is that it is relatively fast when reselecting based on edge type )
        Expanded_edges <- Add_edges_from_mult_sources(edgeData, "info_source", split_string = ";") 
        rm(edgeData)
        edgeList<- color_column2(Expanded_edges, "info_source")
        edgeData <- edgeList[[1]]; edgeLegendData <- edgeList[[2]]
        # unique(edgeData$info_source)
        edgeData <- unique(edgeData)
        # edgeData$rev <- paste0(edgeData$source, edgeData$target)
        #nodeData <-  rbind(sub_net, secondary_nodes) #  "id, name, color, shape"
        network_list <- namedList(AllNodes, edgeData, edgeLegendData)
        network_list
}

remove_free_nodes <- function(rN_ne, nodes, edges){
        if(rN_ne){
                nodesout <- nodes[nodes$id %in% as.vector(unlist(edges[ , c('source', 'target')])) , ]    
        }else{
                nodesout<- nodes
        }   
        nodesout
}
# edges <- edgeData
select_edges <- function(edges, lit = F, met = F, bin = F, reg = F, com = F, kin = F, sig = F){
#        lit <- T; met <- T; bin <- T; reg <- F; com <- F; kin <- T; sig <- F        
       if(lit){lit <- 'literature'}; if(met){met <- 'metablomic'}; if(bin){bin <- "binary"}; if(reg){reg <- 'regulatory'}
       if(com){com <-'complexes'}; if(kin){kin <- 'kinase'}; if(sig){sig <- 'signaling'}
               string <-  paste0(lit,"|" ,met,"|" ,bin,"|" ,reg,"|" ,com, "|" ,kin,"|" ,sig)
      logic <- grepl(string, edges$info_source)
      edges <- edges[logic,]
      edge_key <- data.frame(info_source = c("literature","metabolic","binary","regulatory","complexes","kinase","signaling"),
                                                            color = c("#CC3333","#FF6600","#FFFF33","#00FF00","#00CCCC","#00CCFF","#3366FF"))
      edges <- merge(edges[ , c('source', 'target', 'info_source')], edge_key, by.x = 'info_source')
      edges
}


namedList <- function(...) {
        L <- list(...)
        snm <- sapply(substitute(list(...)),deparse)[-1]
        if (is.null(nm <- names(L))) nm <- snm
        if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
        setNames(L,nm)
}