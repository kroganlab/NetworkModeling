library(plyr)
colA <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n")
colB <- c(1,2,3,1,1,4,4,8,3,1,3,5,6,7)
DF <- data.frame(colA, colB)
list_fa <- unique(DF[duplicated(DF$colB), 'colB' ])
#DF[DF$colB == list_fa[[1]], "colA"]
combination <- combn(DF[DF$colB == list_fa[[1]], "colA"], m = 2)
data1[[1]] <- data.frame(source = combination[1,], target = combination[2,])

data1 <- list()
for(i in 1:length(list_fa)){
  combination <- combn(DF[DF$colB == list_fa[[i]], "colA"], m = 2)
  data1[[i]] <- data.frame(source = combination[1,], target = combination[2,])
}
inter_omics_edges <- rbindlist(data1)

