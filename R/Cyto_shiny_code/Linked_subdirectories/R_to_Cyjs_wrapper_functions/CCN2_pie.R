createCytoscapeNetwork <- 
function(nodeData, edgeData, 
         nodeColor="#888888", nodeShape="ellipse", 
         edgeColor="#888888", edgeSourceShape="none", 
         edgeTargetShape="triangle", nodeHref="") {  
        
        # There must be nodes and nodeData must have at least id and name columns
        if(nrow(nodeData) == 0 || !(all(c("id", "name") %in% names(nodeData)))) {
                return(NULL)
        }
        
        # There must be edges and edgeData must have at least source and target columns
        if(nrow(edgeData) == 0 || !(all(c("source", "target") %in% names(edgeData)))) {
                return(NULL)
        }
        
        # NODES
        ## Add color/shape columns if not present
#         if(!("color" %in% colnames(nodeData))) {
#                 nodeData$color <- rep(nodeColor, nrow(nodeData))
#         }
#         
#         if(!("shape" %in% colnames(nodeData))) {
#                 nodeData$shape <- rep(nodeShape, nrow(nodeData))
#         }
#         
#         if(!("href" %in% colnames(nodeData))) {
#                 nodeData$href <- rep(nodeHref, nrow(nodeData))
#         }
        
        nodeEntries <- NULL
#        nodeData <- nodeData[ , c("id", "name", "nH1N1", "nH3N2", "nH5N1") ]
        for(i in 1:nrow(nodeData)) {   
                tmpEntries <- NULL
                
                for(col in colnames(nodeData)) {
                        if(col == "name" | col == "id" | col == "tim" | col == 'shape' | 
                           col == 'sign_color' | col == "bkgrnd_highlight" ){
                                tmp2 <- paste0(col, ":'", nodeData[i, col], "'")         
                        }else{
                                tmp2 <- paste0(col, ": ", nodeData[i, col], "")    
               #                tmp2 <- paste0(col, ":'", nodeData[i, col], "'")    
                        }
                        tmpEntries <- c(tmpEntries, tmp2)
                }
                
                tmpEntries <- paste(tmpEntries, collapse=", ")
                
                tmp <- paste0("{ data: { ", tmpEntries, "} }")
                
                nodeEntries <- c(nodeEntries, tmp)
        }
        
        nodeEntries <- paste(nodeEntries, collapse=", ")
        
        # EDGES 
        ## Add color/shape columns if not present
        if(!("color" %in% colnames(edgeData))) {
                edgeData$color <- rep(edgeColor, nrow(edgeData))
        }
        
        if(!("sourceShape" %in% colnames(edgeData))) {
                edgeData$edgeSourceShape <- rep(edgeSourceShape, nrow(edgeData))
        }
        
        if(!("targetShape" %in% colnames(edgeData))) {
                edgeData$edgeTargetShape <- rep(edgeTargetShape, nrow(edgeData))
        }
        
        edgeEntries <- NULL
        
        for(i in 1:nrow(edgeData)) {   
                tmpEntries <- NULL
                
                for(col in colnames(edgeData)) {
                        tmp2 <- paste0(col, ":'", edgeData[i, col], "'")
                        tmpEntries <- c(tmpEntries, tmp2)
                }
                
                tmpEntries <- paste(tmpEntries, collapse=", ")
                
                tmp <- paste0("{ data: { ", tmpEntries, "} }")
                
                edgeEntries <- c(edgeEntries, tmp)
        }
        
        edgeEntries <- paste(edgeEntries, collapse=", ")
        
        network <- list(nodes=nodeEntries, edges=edgeEntries)
        
        return(network)
}