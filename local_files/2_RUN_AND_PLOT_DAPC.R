setwd('/Users/chris/Desktop/Workshop_2/Outputs/DAPC/')

#install.packages(adegenet)
library(adegenet)
data<- read.fstat("Lflavomaculatus.dat") 

length(indNames(data)) # number of individuals
grp <- find.clusters(data, max.n.clust=10) # set this for your k range

head(grp$Kstat, 10) # BIC scores per k
grp$stat # Best k based on BIC
head(grp$grp, 59) # Group membership (based on your defined number of clusters)
grp$size # Group size per k cluster

# now do the DAPC
dapc1 <- dapc(data, grp$grp)
myCol <- c("blue","turquoise","lightblue", "darkblue")
scatter(dapc1, scree.da=FALSE, scree.pca=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=1, cex=1.5, clab=0, leg=FALSE, txt.leg=paste("Cluster",1:4))
png("../DAPC/plot_Lflavomaculatus_DAPC_k=4.png",width=500,height=500)
dev.off()


