
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~ K-MEANS ~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit = kmeans(dat.wide, centers=cluster_size, iter.max = 1000, nstart=100)

x.kmeans = cbind(dat.wide, cluster=fit$cluster)
x.kmeans$id <- row.names(x.kmeans)
x.kmeans <- melt(x.kmeans, id=c("id",'cluster'), variable.name='time', value.name='log2fc')



p <- ggplot(x.kmeans, aes(x=time, y=log2fc, group=id)) + geom_line(aes(color=cluster)) + facet_wrap(~cluster, ncol=5)
ggsave('~/Desktop/kmeans.pdf', p, height=20, width=10)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~ K-NN ~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


