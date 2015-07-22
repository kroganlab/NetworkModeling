setwd("~/Box Sync/Doug/projects/FLU_Networking_code/Cyto_shiny_code")
library(shiny)
library(shinyapps)

shinyapps::configureApp("Cyto_shiny_code", size="large")

terminateApp("CytoShiny_flu")

source("http://bioconductor.org/biocLite.R")
biocLite("DLBCL")

library(BioNet)
library(DLBCL)
data(dataLym)
data(interactome)

?interactome

pvals <- cbind(t = dataLym$t.pval, s = dataLym$s.pval)

rownames(pvals) <- dataLym$label
pval <- aggrPvals(pvals, order = 2, plot = FALSE)

?subNetwork

subnet <- subNetwork(dataLym$label, interactome)
subnet <- rmSelfLoops(subnet)
?rmSelfLoops
subnet

fb <- fitBumModel(pval, plot = FALSE)
scores <- scoreNodes(subnet, fb, fdr = 0.001)

?scoreNodes

module <- runFastHeinz(subnet, scores)
logFC <- dataLym$diff
names(logFC) <- dataLym$label

plotModule(module, scores = scores, diff.expr = logFC)


library(igraph)
el <- cbind(c("a", "b", "c", "d", "e", "f", "d"), c("b", "c", "d", "e", "f", "a", "b"))
graph <- graph.edgelist(el, directed=TRUE)

?graph.edgelist

node.list <- c("a", "b", "c")
graph2 <- subNetwork(nodeList=node.list, network=graph)
## Not run: par(mfrow=c(1,2));
plotModule(graph);
plotModule(graph2)

# or in graphNEL format: 
graph3 <- igraph.to.graphNEL(graph)
graph4 <- subNetwork(nodeList=node.list, network=graph3)
## End(Not run)