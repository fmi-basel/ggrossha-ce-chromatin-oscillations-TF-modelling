# This R script contains code that reproduces the timing tranjectory plots for the wildtype RNA-seq time course experiment

# read gene expression count table (GSE288875 on GEO:https://www.ncbi.nlm.nih.gov/geo/)
CO <- read.delim("data/GSE288875_expr_mRNA_rawCounts.txt",check.names=FALSE)
CON <- t(t(CO[,-1])/colSums(CO[,-1])*1000000) # remove first column (gene size) and normalize to 1M for timing plot

# split the sample names containing annotation information. this needs to be adjusted for every dataset
annot <- do.call(rbind,strsplit(colnames(CON),"_"))

# split the dataset according to the various conditions. this needs to be adjusted for every dataset
condL <- split(1:nrow(annot),annot[,3])

TL <- log2(CON+1) # log2 transform after adding a pseudo count of 1

# mean normalization for each gene. this needs to be adjusted for every dataset. in this case, simply the mean of all the
# samples is subtracted. in the case of treatment/mutant conditions, it is best to substract the mean from one or more complete
# cycles in the non-perturbed condition. failure of proper mean normalization will likely produce timing trajectories that
# that are not centered around (x,y)=(0,0)
TLD <- TL - rowMeans(TL)

# load a table with oscillating genes and associated phases, synchronized to larval cycle (MoltSync). 0 degree refers to the
# beginning of a larval stage and 280 degrees refers roughly to the start of the molt (GSE288873 on GEO:https://www.ncbi.nlm.nih.gov/geo/)
oscGenes <-  read.delim("data/GSE288873_AllOsc_info_CIclass_MoltSync.txt")

# split the genes according to their phase bin (0,90,180,270)
oscGenesBinGenesL <- split(rownames(oscGenes),oscGenes$GeneSet)

# caluculate average expression for each phase bin
oscGenesBinExpr <- sapply(oscGenesBinGenesL,function(x){colMeans(TLD[rownames(TLD) %in% x,])})

# calculate the two orthogonal axes
oscGenesPR <- cbind(X=oscGenesBinExpr[,"GeneSet_0"]-oscGenesBinExpr[,"GeneSet_180"],Y=oscGenesBinExpr[,"GeneSet_90"]-oscGenesBinExpr[,"GeneSet_270"])

#define color palette for the samples
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# plot the timing trajectories for each condition separately
#pdf("plots/timing_MoltSync.pdf",width=18,height=8,pointsize=20)
par(mar=c(4,4,2,4),mfrow=c(1,length(condL)),xpd=TRUE)
plRan <- c(-max(abs(oscGenesPR)),max(abs(oscGenesPR)))
for(i in 1:length(condL)){
  plot(oscGenesPR[condL[[i]],],xlim=plRan,ylim=plRan,col="darkgray",type='l',main=names(condL)[i])
  points(oscGenesPR[condL[[i]],],pch='.',cex=18,col=jet.colors(length(condL[[i]])))
  legend("topright",inset=c(-0.21,0),legend=annot[condL[[1]],2],fill=jet.colors(length(condL[[i]])),bty="n")
}
#dev.off()

