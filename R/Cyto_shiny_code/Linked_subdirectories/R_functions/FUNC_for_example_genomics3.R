Add_edges_from_mult_sources <- function(edges, source_col, split_string = ";"){
        # edges <- edgeData; source_col <- "info_source"; split_string <- ';'
        k <- 1
        add_edges <- edges[FALSE, ]
        for(i in 1:nrow(edges)){
     #           i <- 7
                sources <- unlist(strsplit(gsub(split_string,"_" ,edges[i, source_col]), as.vector("_")))
                num_sources <- length(sources)
                if(num_sources > 1){
                        for(j in 1:num_sources){
                                add_edges[k,] <- edges[i, ]
                                add_edges[k, source_col ] <- sources[j]
                                k <- k+1
                        }
                }
        }
        cleaned_edges <- edges[!grepl(split_string, edges[ , source_col]),]
        new_edges <- rbind(cleaned_edges, add_edges)
}



# wrapper functions for RCy3 package functions
get_node_strains_attributes <- function(M){
        M$strains <- paste0(ifelse(M$criteria.H5N1 == 1,"H5N1Z",""),
                            ifelse(M$criteria.H3N2 == 1,"ZH3N2Z",""),
                            ifelse(M$criteria.H1N1 == 1,"ZH1N1",""))
        M$strains <- gsub("ZZ", "_", M$strains) 
        M$strains <- gsub("Z", "", M$strains)
        M
}

set_node_attributes_Cyto <- function(M,colname,g_in, type){
        attribute_values <- unique(M[ , colname])
        for(i in 1:length(attribute_values)){
                nodeData (g_in, as.character(M[M[,colname] == attribute_values[i],
                                               type] ), colname ) <- attribute_values[i]
        }
        g_in
}
##
color_column3 <- function(node_data, column_name){
#         node.colors.18 <- c("#CC3333", "#FF6600", "#3366FF", "#00FF00", "#00CCCC","#00CCFF" ,"#669955",
#                             "#9933FF", "#FF00FF", "#FFCCCC", "#FFCC99", "#FFFFCC", "#CCFFCC", "#99FFCC",
#                             "#CCFFFF" ,"#99CCFF", "#CCCCFF", "#FFCCFF" )
        node.colors.18 <- c("#FF0000", "#000000", "#0000FF", "008000", "#800000","#8D38C9" ,"#F88017")
#                             "#9933FF", "#FF00FF", "#FFCCCC", "#FFCC99", "#FFFFCC", "#CCFFCC", "#99FFCC",
#                             "#CCFFFF" ,"#99CCFF", "#CCCCFF", "#FFCCFF" )
        tags_for_colors <- as.character(unique(node_data[ , column_name]))
        node.colors <- node.colors.18[1:length(tags_for_colors)]
        colors_and_tags <- data.frame(cbind(node.colors, tags_for_colors))
        colnames(colors_and_tags) <- c('color', 'info_source')
        node_data2 <- merge(node_data, colors_and_tags, by='info_source')   
        colnames(colors_and_tags) <- c('node.colors', 'tags_for_colors')
        list(node_data2, colors_and_tags)
}
# color_column2 <- function(node_data, column_name){
#          node.colors.18 <- c("#CC3333", "#FF6600", "#3366FF", "#00FF00", "#00CCCC","#00CCFF" ,"#669955",
#                              "#9933FF", "#FF00FF", "#FFCCCC", "#FFCC99", "#FFFFCC", "#CCFFCC", "#99FFCC",
#                              "#CCFFFF" ,"#99CCFF", "#CCCCFF", "#FFCCFF" )
#          node.colors.18 <- c("#FF0000", "#000000", "#0000FF", "008000", "#800000","#8D38C9" ,"#F88017",
#                              "#9933FF", "#FF00FF", "#FFCCCC", "#FFCC99", "#FFFFCC", "#CCFFCC", "#99FFCC",
#                              "#CCFFFF" ,"#99CCFF", "#CCCCFF", "#FFCCFF" )
#         tags_for_colors <- as.character(unique(node_data[ , column_name]))
#         node.colors <- node.colors.18[1:length(tags_for_colors)]
#         colors_and_tags <- data.frame(cbind(node.colors, tags_for_colors))
#         colnames(colors_and_tags) <- c('node.colors', 'tags_for_colors')
#         for(i in 1:nrow(node_data)){
#                 node_data[i , "color"] <- as.character(colors_and_tags[colors_and_tags$tags_for_colors == 
#                                                                                as.character(node_data[i, column_name]) , "node.colors"])  
#         }
#         list(node_data, colors_and_tags)
# }
##
##



?data.frame
color_column <- function(node_data, column_name){
        node.colors.18 <- c("#CC3333", "#FF6600", "#FFFF33", "#00FF00", "#00CCCC","#00CCFF" ,"#3366FF",
                            "#9933FF", "#FF00FF", "#FFCCCC", "#FFCC99", "#FFFFCC", "#CCFFCC", "#99FFCC",
                            "#CCFFFF" ,"#99CCFF", "#CCCCFF", "#FFCCFF" )
        tags_for_colors <- unique(node_data[ , column_name])
        node.colors <- node.colors.18[1:length(tags_for_colors)]
        colors_and_tags <- data.frame(cbind(node.colors, tags_for_colors))
        for(i in 1:nrow(node_data)){
                node_data[i , "color"] <- as.character(colors_and_tags[colors_and_tags$tags_for_colors == 
                                                                               node_data[i, column_name] , "node.colors"])      
                
        }
        node_data
}


remove_dup_edge_rev_order <- function(edges_sort, colname1, colname2, attribute_colname1){
# edges_sort <- edges; colname1 <- 'source'; colname2 <- 'target'; attribute_colname1 <- 'info_source'
        for(i in 1:nrow(edges_sort)){
                if(edges_sort[i, colname2] < edges_sort[i,colname1]){
                        edges_sort[i, "target2"] <- edges_sort[i, colname1]
                        edges_sort[i, "source2"] <- edges_sort[i, colname2] 
                } else{
                        edges_sort[i, "target2"] <- edges_sort[i, colname2]
                        edges_sort[i, "source2"] <- edges_sort[i, colname1] 
                }
        }
        edges_sort2 <- unique(edges_sort[, c("source2","target2", attribute_colname1)])
        
        colnames(edges_sort2) <- c("source", "target", attribute_colname1)
        edges_sort2
} # end function definition

Prey_connections <- function(Assoc_DF, col1, col2, List_Prey, col3){
        Assoc_list <- Assoc_DF[Assoc_DF[ , col1]  %in% List_Prey[ , col3] & 
                                       Assoc_DF[ , col2] %in% List_Prey[ , col3] , c(col1, col2)]
        Assoc_list <- unique(Assoc_list)
}
Prey_connections_OR3 <- function(Assoc_DF, col1, col2, List_Prey, col3, sourcecol){
        Assoc_list1 <- Assoc_DF[Assoc_DF[ , col1]  %in% List_Prey[ , col3]  , c(col1, col2, sourcecol)]
        Assoc_list2 <- Assoc_DF[Assoc_DF[ , col2]  %in% List_Prey[ , col3]  , c(col2, col1, sourcecol)]
        colnames(Assoc_list2) <- c(col1, col2 , sourcecol)
        Assoc_list <- rbind(Assoc_list1,Assoc_list2)
        Assoc_list <- unique(Assoc_list)
}

Prey_connections_OR <- function(Assoc_DF, col1, col2, List_Prey, col3){
        Assoc_list <- Assoc_DF[Assoc_DF[ , col1]  %in% List_Prey[ , col3] & 
                                       Assoc_DF[ , col2] %in% List_Prey[ , col3] , c(col1, col2)]
        Assoc_list <- unique(Assoc_list)
}

Prey_connections_OR2 <- function(Assoc_DF, col1, col2, List_Prey, col3){
        Assoc_list1 <- Assoc_DF[Assoc_DF[ , col1]  %in% List_Prey[ , col3]  , c(col1, col2)]
        Assoc_list2 <- Assoc_DF[Assoc_DF[ , col2]  %in% List_Prey[ , col3]  , c(col2, col1)]
        colnames(Assoc_list2) <- c(col1,col2)
        Assoc_list <- rbind(Assoc_list1,Assoc_list2)
        Assoc_list <- unique(Assoc_list)
}



Prey_connections_2vec <- function(Assoc_DF, col1, col2, List_Prey1, List_Prey2, col3_1, col3_2){
        Assoc_list1 <- Assoc_DF[Assoc_DF[ , col1]  %in% List_Prey1[ , col3_1] & 
                                       Assoc_DF[ , col2] %in% List_Prey2[ , col3_2] , c(col1, col2)]
        Assoc_list2 <- Assoc_DF[Assoc_DF[ , col2]  %in% List_Prey1[ , col3_1] & 
                                        Assoc_DF[ , col1] %in% List_Prey2[ , col3_2] , c(col1, col2)]
        Assoc_list <- rbind(Assoc_list1,Assoc_list2)
        Assoc_list <- unique(Assoc_list)
}

Produce_Venn <- function(M, time, idtype ){
        M <- unique(M[M$time == time, c("strain", idtype, "criteria") ])
        DF_reshape <- reshape(M, timevar = "strain" , 
                              idvar = idtype, direction = "wide")
        DF_reshape[is.na(DF_reshape)] <- 0
        names = gsub("criteria.", "", colnames(DF_reshape[,c(2:ncol(DF_reshape))]))
        vennCountTable <- vennCounts(DF_reshape[,c(2:ncol(DF_reshape))])
        vennPackage<- list(vennCountTable = vennCountTable, names = names)
        vennPackage
}
sub_get <- function(M, time, idtype ){
#        M <- unique(M[M$time == time, c("strain", idtype, "criteria", "symbol") ])
        M <- unique(M[M$time == time, c("strain", idtype, "criteria") ])
        DF_reshape <- reshape(M, timevar = "strain" , 
                              idvar = idtype, direction = "wide")
        DF_reshape[is.na(DF_reshape)] <- 0
        names = gsub("criteria.", "", colnames(DF_reshape[,c(2:ncol(DF_reshape))]))
        DF_reshape 
}

sub_get_strain <- function(M, STRAIN, idtype ){
        M <- unique(M[M$strain == STRAIN, c("time", idtype, "criteria") ])
        DF_reshape <- reshape(M, timevar = "time" , 
                              idvar = idtype, direction = "wide")
        DF_reshape[is.na(DF_reshape)] <- 0
        names = gsub("criteria.", "", colnames(DF_reshape[,c(2:ncol(DF_reshape))]))
        DF_reshape 
}

Produce_Venn2 <- function(M, STRAIN , idtype ){
        #M <- gm_sub; STRAIN <- "H1N1" ; idtype <- "uniprot_ac"
        M <- unique(M[M$strain == STRAIN, c("time", idtype, "criteria") ])
        DF_reshape <- reshape(M, timevar = "time" , 
                              idvar = idtype, direction = "wide")
        DF_reshape[is.na(DF_reshape)] <- 0
        names = gsub("criteria.", "", colnames(DF_reshape[,c(2:ncol(DF_reshape))]))
        vennCountTable <- vennCounts(DF_reshape[,c(2:ncol(DF_reshape))])
        vennPackage <- list(vennCountTable = vennCountTable, names = names)
        vennPackage
}


plotVenn <- function(M, time, idtype){
        vennObject <- Produce_Venn(M, time, idtype)
        vennDiagram(vennObject$vennCountTable, names = vennObject$names, 
                    main = paste0(time, " ", idtype))
}

plotVenn_times <- function(M, STRAIN, idtype){
        vennObject <- Produce_Venn2(M, STRAIN, idtype)
        #par(mfrow = c(1,1))
        vennDiagram(vennObject$vennCountTable, names = vennObject$names,
                    main = paste0(STRAIN, " ", idtype))
}
Get_Names2 <- function(M, time, Strains, idtype){
        OUT <- 0
        M <- unique(M[M$time == time, c("strain", idtype, "criteria") ])
        DF_reshape <- reshape(M, timevar = "strain" , 
                              idvar = idtype, direction = "wide")
        DF_reshape[is.na(DF_reshape)] <- 0
        colnames(DF_reshape) = c(idtype, gsub("criteria.", "", colnames(DF_reshape[,c(2:ncol(DF_reshape))])))
        if(length(Strains) == 3){OUT <- DF_reshape[((DF_reshape[, Strains[1]] + DF_reshape[, Strains[2]] + 
                                                             DF_reshape[,Strains[3]]) == length(Strains)), idtype]}
        if(length(Strains) == 2){OUT <- DF_reshape[((DF_reshape[, Strains[1]] + 
                                                             DF_reshape[, Strains[2]])  == length(Strains)), idtype]}
        OUT
}

Get_Names2F <- function(M, time, Strains, idtype){
        OUT <- 0
        M <- unique(M[M$time == time, c("strain", idtype, "criteria", "log2fc", "q_value") ])
        DF_reshape <- reshape(M, timevar = "strain" , 
                              idvar = idtype, direction = "wide")
        DF_reshape[is.na(DF_reshape)] <- 0
        colnames(DF_reshape) = c(idtype, gsub("criteria.", "", colnames(DF_reshape[,c(2:ncol(DF_reshape))])))
        if(length(Strains) == 3){OUT <- DF_reshape[((DF_reshape[, Strains[1]] + DF_reshape[, Strains[2]] + 
                                                             DF_reshape[,Strains[3]]) == length(Strains)), idtype]}
        if(length(Strains) == 2){OUT <- DF_reshape[((DF_reshape[, Strains[1]] + 
                                                             DF_reshape[, Strains[2]])  == length(Strains)), idtype]}
        OUT
}
assembleDF <- function(M, idtype){
        DF <- data.frame(); n <- 1
        times <- c("12H","01D","02D","03D","04D")
        Strains <- list(c("H1N1","H5N1", "H3N2") ,c("H1N1", "H5N1") , c("H1N1" ,"H3N2" ) ,c("H5N1", "H3N2") )
        Strains_set_names <- c("H1N1_H5N1_H3N2","H1N1_H5N1","H1N1_H3N2","H5N1_H3N2")
        for(i in 1:4){ # Strain combination loop
                for(j in 1:5){ # time loop
                        entry_list <- Get_Names2(M, times[j], Strains[[i]], idtype)
                        for(k in 1:length(entry_list)){ # fill in rows loop
                                DF[n , "time"] <- times[j]
                                DF[n , "strains"] <- Strains_set_names[i]
                                DF[n , idtype] <- entry_list[k]
                                n <- n+1
                        }
                }
        }
        DF
}


# Assoc_DF <- Mouse_Map; col1 <- "geneid_1"; col2 <- "geneid_1"
#         Assoc_list1 <- Assoc_DF[Assoc_DF[ , col1]  %in% List_Prey1[ , col3_1] & 
#                                         Assoc_DF[ , col2] %in% List_Prey2[ , col3_2] , c(col1, col2)]
#         Assoc_list2 <- Assoc_DF[Assoc_DF[ , col2]  %in% List_Prey1[ , col3_1] & 
#                                         Assoc_DF[ , col1] %in% List_Prey2[ , col3_2] , c(col1, col2)]
#         Assoc_list <- rbind(Assoc_list1,Assoc_list2)
#         Assoc_list <- unique(Assoc_list)
# }