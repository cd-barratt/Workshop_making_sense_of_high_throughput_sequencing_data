setwd('/Users/chris/Desktop/RAD-seq_workshop/Workshop_2/Outputs/Admixture/')

# plot CV scores
results <- read.delim('/Users/chris/Desktop/RAD-seq_workshop/Workshop_2/work/Lflavomaculatus/Outputs/Admixture/RESULTS_Lflavomaculatus_admixture_cleaned.txt', header=FALSE)
plot (results, pch=19, col= "black", type = "b", cex=0.75, lty=1, xlab="k clusters", ylab="CV score", main = "Cross validation scores across k")
dev.copy(tiff,"./RESULTS_CV_scores")
dev.off()

# plot bars
tbl=read.table("/Users/chris/Desktop/Workshop_2/Outputs/Admixture/Lflav.4.Q")
ord = tbl[order(tbl$V1,tbl$V3,tbl$V2,tbl$V4),]

labels <-read.delim("../../Input_files/Admixture/popinfo_labels.txt", header=FALSE)
labels <- data.frame(labels[,2])
ord_with_labels <- data.frame(ord[,1:4], row.names=labels[,1])
par(mar=c(12,4,1,1))
bp = barplot(t(as.matrix(ord_with_labels)), space = (0.04), col=c("blue", "turquoise", "darkblue", "lightblue"), xlab="", ylab="Ancestry", border=1, las=2, cex.names=0.9)
dev.copy(pdf,"simple_plot_labelled_individuals.pdf", width=30, height=6)
dev.off()

labels <-read.delim("../../Input_files/Admixture/popinfo_labels.txt", header=FALSE)
labels <- data.frame(labels[,5])
ord_with_labels <- data.frame(ord[,1:4], row.names=labels[,1])
par(mar=c(12,4,1,1))
bp = barplot(t(as.matrix(ord_with_labels)), space = (0.04), col=c("blue", "turquoise", "darkblue", "lightblue"), xlab="", ylab="Ancestry", border=1, las=2, cex.names=0.9)
dev.copy(pdf,"simple_plot_labelled_localities.pdf", width=30, height=6)
dev.off()
