# This R script contains code that reproduces the various linear models used to predict ATAC-seq amplitude and phase from Transcrition factor (TF) ChIP data

# plot circles with radii x
plotCircles <- function(x,...){
  phivec <- seq(0,2*pi,length.out=100)
  for(i in 1:length(x)){lines(x[i]*cbind(cos(phivec),sin(phivec)),type='l',...)}
}


# plot radar chart with circles with radii x
plotRadarAnnot <- function(x,...){
  phivec <- seq(0,2*pi,length.out=100)
  phivecs <- 2*pi/8*c(0:7)
  radSegs <- max(x)*exp(1i*phivecs)
  for(i in 1:length(x)){lines(x[i]*cbind(cos(phivec),sin(phivec)),...)}
  segments(0,0,Re(radSegs),Im(radSegs),...)
  text(1.1*Re(radSegs),1.1*Im(radSegs),paste0(phivecs/pi*180,"\u00B0"))
}

# load master table with all the required information (GSE288867 on GEO:https://www.ncbi.nlm.nih.gov/geo/)
MT <- read.delim("data/GSE288867_ATAC_peaks_xy_linModelCont_allInfo.txt",check.names=FALSE)

# extract the ChIP-seq enrichment columns
ChIP <- as.matrix(MT[,grep("ChIP_",colnames(MT))])
colnames(ChIP) <- gsub("ChIP_","",colnames(ChIP))

# extract the coefficients from the cosine fit performed on the ATAC-seq peaks (x=ampl*cos(phi),y=ampl*sin(phi))
ATAC_xy <- as.matrix(MT[,c("ATAC_x","ATAC_y")])

# plot the results from the cosine fit (x=ampl*cos(phi),y=ampl*sin(phi)) performed on the ATAC-seq data. included are all oscillating 
# ATAC-seq peaks as well as control regions that do not oscillate but are only accessible only in oscillating tissues (see paper)
#pdf("plots/predict_ATAC_by_currTFs_xy.pdf",pointsize=15)
plot(ATAC_xy,col=c("black","darkgray")[(MT$class != "osc") + 1],xlab="x",ylab="y",pch='.',cex=3,xlim=c(-1.7,1.7),ylim=c(-1.7,1.7))
plotRadarAnnot(c(0.5,1,1.5),col="darkgray")
legend('bottomleft', legend = c("Oscillating","Not oscillating, only accessible in oscillating cell types"), pch = 16, col = c("black","darkgray"),bty="n")
#dev.off()

# fit two linear models, one that predicts x and one that predicts y from ChIP-seq enrichments, no intercept is used
fm_x <- lm(ATAC_xy[,1] ~ ChIP + 0)
fm_y <- lm(ATAC_xy[,2] ~ ChIP + 0)

# extract the coefficients from the model into a vector with complex numbers
fm_coeffs_cplx <- fm_x$coefficients + 1i*fm_y$coefficients
names(fm_coeffs_cplx) <- colnames(ChIP)

# for every ATAC-seq peak, calculate the contribution vectors for all TFs
fmPredM <- t(t(ChIP) * fm_coeffs_cplx)

# define color spaces
vir.colors <- colorRampPalette(c("#440154FF","#46337EFF","#365C8DFF","#277F8EFF","#1FA187FF","#4AC16DFF","#9FDA3AFF","#FDE725FF"))
rnb.colors <- colorRampPalette(c("#FF5400","#FF8E00","#FFD200","#81E650","#00D267","#00C0FF","#8B48FE","#CA41FC","#FF46FB"))

# define a color for each TF
tf_cols <- rnb.colors(ncol(ChIP))
names(tf_cols) <- colnames(ChIP)

lm_adjRsqLab <-paste0(c("x:","y:"),round(c(summary(fm_x)$adj.r.squared,summary(fm_y)$adj.r.squared),3))

# plot the coefficients from the two linear models
#pdf("plots/predict_ATAC_by_currTFs_coeffs_lin.pdf",pointsize=16,width=7,height=7.7)
plot(fm_x$coefficients,fm_y$coefficients,col=tf_cols,pch=1,cex=2,lwd=2.5,xlim=c(-0.29,0.29),ylim=c(-0.29,0.29),xlab="Model coefficient x\'",ylab="Model coefficient y\'")
segments(0,0,fm_x$coefficients,fm_y$coefficients,lwd=2.5,col=tf_cols)
text(fm_x$coefficients,fm_y$coefficients,colnames(ChIP),pos=2)
plotCircles(c(0.1,0.2,0.3),lty=3,col="darkgray")
legend('topleft', legend = c(lm_adjRsqLab),title = "Adj R^2")
#dev.off()

# bin the ChIP-seq enrichment values
PP_cat_breaks <- c(-100,seq(-0.25,2.25,by=0.5),100)
PP_cat <- apply(ChIP,2,function(x){as.character(as.numeric(cut(x,PP_cat_breaks)))})

# fit two linear models using binned ChIP-seq enrichment data as predictors
fm_x_cat <- lm(ATAC_xy[,1] ~  . ,data=data.frame(PP_cat))
fm_y_cat <- lm(ATAC_xy[,2] ~  . ,data=data.frame(PP_cat))

fm_cat_dat <- data.frame(TF=colnames(ChIP)[fm_y_cat$assign[-1]],x=fm_x_cat$coefficients[-1],y=fm_y_cat$coefficients[-1])
fm_cat_dat_L <- split(fm_cat_dat,fm_cat_dat[,1])

lm_cat_adjRsqLab <-paste0(c("x:","y:"),round(c(summary(fm_x_cat)$adj.r.squared,summary(fm_y_cat)$adj.r.squared),3))

pchDia <- (1:(nrow(fm_cat_dat_L[[1]])+1))^2/24+11/24

#pdf("plots/predict_ATAC_by_currTFs_coeffs_linLevels.pdf",pointsize=16,width=7,height=7.7)
par(xaxs="i",yaxs="i")
plot(rep(0.75,length(pchDia)),seq(-0.75,-0.15,length.out=length(pchDia)),xlab="x",ylab="y",xlim=c(-0.8,0.8),ylim=c(-0.8,0.8),lwd=2,cex=pchDia)
text(rep(0.67,length(pchDia)),seq(-0.75,-0.15,length.out=length(pchDia)),1:length(pchDia))
for(i in 1:length(fm_cat_dat_L)){
  lines(rbind(c(0,0),fm_cat_dat_L[[colnames(ChIP)[i]]][,2:3]),type='o',lwd=2,col=tf_cols[i],pch=1,cex=pchDia)
}
plotCircles(c(0.2,0.4,0.6,0.8),lty=3,col="darkgray")
fLab <- do.call(rbind,lapply(fm_cat_dat_L,function(x){x[nrow(x),2:3]}))
text(fLab[,1]-0.025,fLab[,2],rownames(fLab),pos=2)
legend('topleft', legend = c(lm_cat_adjRsqLab),title = "Adj R^2")
#dev.off()


# add interaction terms to the two linear models
fm_x_int <- lm(ATAC_xy[,1] ~ 0 + .^2, data=data.frame(ChIP,check.names=FALSE))
fm_y_int <- lm(ATAC_xy[,2] ~ 0 + .^2, data=data.frame(ChIP,check.names=FALSE))

# extract a combined p-value for every predictor considering both models
fm_x_int_p <- summary(fm_x_int)$coefficients[,"Pr(>|t|)"]
fm_y_int_p <- summary(fm_y_int)$coefficients[,"Pr(>|t|)"]
fm_int_comb_p <- apply(cbind(fm_x_int_p,fm_y_int_p),1,function(p){pchisq((sum(log(p))*-2), df=length(p)*2, lower.tail=F)})
names(fm_int_comb_p) <- gsub("`","",names(fm_int_comb_p))


lm_int_adjRsqLab <-paste0(c("x:","y:"),round(c(summary(fm_x_int)$adj.r.squared,summary(fm_y_int)$adj.r.squared),3))

cIntCol <- c("black","red")[c(rep(1,ncol(ChIP)),rep(2,length(fm_int_comb_p)-ncol(ChIP)))]

# plot the p values for all coefficients from the model which includes interaction terms
#pdf("plots/predict_ATAC_by_currTFs_coeffs_Interactions.pdf",pointsize=11,width=8,height=4)
par(mar=c(8,4,2,2))
barplot(-log10(fm_int_comb_p),las=2,ylab="combined p-value (-log10)",col=cIntCol)
legend(x="topright",legend=c("Individual linear predictors","Interaction terms"),fill=c("black","red"),bty="n")
legend('right', legend = c(lm_int_adjRsqLab),title = "Adj R^2")
#dev.off()


# define the training and test set
selTR <- MT$training
selTE <- !selTR

# fit models using only the test data
fm_x_tr <- lm(ATAC_xy[selTR,1] ~ ChIP[selTR,] + 0)
fm_y_tr <- lm(ATAC_xy[selTR,2] ~ ChIP[selTR,] + 0)
fm_coeffs_tr_cplx <- fm_x_tr$coefficients + 1i*fm_y_tr$coefficients
names(fm_coeffs_tr_cplx) <- colnames(ChIP)

fmPredM_tr <- t(t(ChIP) * fm_coeffs_tr_cplx)

pLabT <- (rowMeans(cbind(fm_x$coefficients,fm_x_tr$coefficients))+1i*rowMeans(cbind(fm_y$coefficients,fm_y_tr$coefficients))) * exp(1i*-0.0)

# plot the results from the model that uses all the data and the one that uses only the training data
#pdf("plots/predict_ATAC_by_currTFs_coeffs_lin_All_vs_Train10.pdf",pointsize=16,width=7,height=7.7)
plot(fm_x$coefficients,fm_y$coefficients,col=tf_cols,pch=1,cex=2,lwd=2,xlim=c(-0.29,0.29),ylim=c(-0.29,0.29),xlab="Model coefficient x",ylab="Model coefficient y")
segments(0,0,fm_x$coefficients,fm_y$coefficients,lwd=2.5,col=tf_cols)
points(fm_x_tr$coefficients,fm_y_tr$coefficients,col=tf_cols,pch=16,cex=2,lwd=2)
text(Re(pLabT)-0.01,Im(pLabT)-0.01,colnames(ChIP),pos=2)
plotCircles(c(0.1,0.2,0.3),lty=3,col="darkgray")
abline(h=0,lty=3,col="darkgray")
abline(v=0,lty=3,col="darkgray")
legend('topleft', legend = c("All the data","10% of the data"), pch = c(1,16),bty="n")
#dev.off()


# fetch TF cluster information
km_res_clu <- setNames(MT$TFCluster,rownames(MT))

# calculate average TF enrichments for all the clusters (only the test data is used here, not the traing data!)
ChIP_cl <- aggregate(ChIP[selTE,],list(km_res_clu[selTE]),mean)
rownames(ChIP_cl) <- ChIP_cl[,1]
ChIP_cl <- as.matrix(ChIP_cl[,-1])

# determine variance explained in ChIP enrichment signal by the clusters
ChIP_clPred <- sapply(1:ncol(ChIP),function(i){lm(ChIP[,i] ~ as.character(km_res_clu))$fitted.values})
ChIP_var_explained_by_clusters <- cor(as.vector(ChIP),as.vector(ChIP_clPred))^2

# split the ChIP-seq enrichment matrix into multiple panels for visual clarity
panelToCluL <- split(rownames(ChIP_cl),0:(nrow(ChIP_cl)-1) %/% round(nrow(ChIP_cl)/5))

#pdf("plots/predict_ATAC_by_currTFs_TfClusters.pdf",pointsize=14)
pZrange <- c(-0.5,3.5)
par(mfrow=c(1,length(panelToCluL)),mar=c(6,0.5,0.5,2),xpd=TRUE)
for(i in 1:length(panelToCluL)){
  peaksInPanelM <- do.call(rbind,lapply(as.numeric(panelToCluL[[i]]),function(x){cbind(x,rownames(ChIP)[km_res_clu==x])}))
  tM <- pmin(pmax(ChIP[peaksInPanelM[,2],],pZrange[1]),pZrange[2])
  rLab <- t(sapply(split(seq(1,0,length.out=nrow(tM)),peaksInPanelM[,1]),range))
  cLab <- rowMeans(rLab)

  image(t(tM[nrow(tM):1,]),col=vir.colors(256),axes=FALSE,useRaster=TRUE,zlim=pZrange)
  box()
  axis(1,seq(0,1,length.out=ncol(tM)),colnames(tM),las=2)
  text(1.22,cLab,names(cLab),las=2)
  segments(1.12,rLab[,1]+0.003,1.12,rLab[,2]-0.003)
  segments(1.09,rLab[,1],1.12,rLab[,1]+0.003)
  segments(1.09,rLab[,2],1.12,rLab[,2]-0.003)
}
#dev.off()


# create a separate plot for the colorbar
#pdf("plots/predict_ATAC_by_currTFs_TfClusters_colorbar.pdf",width=4,height=2.3)
image(seq(pZrange[1],pZrange[2],length=256),1,matrix(1:256),col=vir.colors(256),yaxt='n',xlab="Normalized ChIP enrichment (log2)",ylab="")
#dev.off()

# calculate the contribution vectors for all TFs at the TF cluster level
ChIP_cl_comp <- t(t(ChIP_cl) * fm_coeffs_tr_cplx)

# calculate the average x and y at for each TF cluster (only using the test data!)
MB_real <- aggregate(ATAC_xy[selTE,],list(km_res_clu[selTE]),mean)
rownames(MB_real) <- MB_real[,1]
MB_real <- as.matrix(MB_real[,-1])

# calculate the predicted average x and y at for each TF cluster (only using the test data!)
MB_pred <- aggregate(cbind(rowSums(Re(fmPredM_tr)),rowSums(Im(fmPredM_tr)))[selTE,],list(km_res_clu[selTE]),mean)
rownames(MB_pred) <- MB_pred[,1]
MB_pred <- as.matrix(MB_pred[,-1])

# compare the measured x and y for each TF cluster to the predicted
#pdf("plots/predict_ATAC_by_currTFs_additive_TfClusters_allPairsScatter.pdf")
plot(MB_real,pch=16,xlab="x",ylab="y",xlim=c(-0.7,0.7),ylim=c(-0.7,0.7))
points(MB_pred,col="red",pch=16)
segments(MB_real[,1],MB_real[,2],MB_pred[,1],MB_pred[,2],col="darkgray",lwd=2)
plotCircles(c(0.2,0.4,0.6),lty=3,col="darkgray")
legend('topright', legend = c("Measured","Predicted"), pch = 16, col = 1:2, ,bty="n")
legend('bottomleft', legend = c(paste0("x:",round(cor(MB_real[,1],MB_pred[,1]),3)),paste0("y:",round(cor(MB_real[,2],MB_pred[,2]),3))),title = "r")
#dev.off()

# select the clusters to plot in main figure
clusToPlot1 <- sort(order(rowSums(MB_pred^2),decreasing=TRUE)[1:24])

#pdf("plots/predict_ATAC_by_currTFs_additive_TfClusters.pdf",pointsize=18,width=16,height=13)
par(mfrow=c(5,6),mar=c(0.5,0.5,0.5,0.5))
for(i in clusToPlot1){
  plot(0,0,xlim=c(-0.7,0.7),ylim=c(-0.7,0.7),type="n",xlab="",ylab="",axes=FALSE)
  axis(1,labels = FALSE)
  axis(2,labels = FALSE)
  box()
  plotCircles(c(0.2,0.4,0.6),lty=2,col="gray")

  fOrd <- order(abs(ChIP_cl_comp[i,]))
  seM <- cbind(c(0,ChIP_cl_comp[i,fOrd][-length(fOrd)]),ChIP_cl_comp[i,fOrd])
  tCols <- tf_cols[fOrd]
  segments(cumsum(Re(seM[,1])),cumsum(Im(seM[,1])),cumsum(Re(seM[,2])),cumsum(Im(seM[,2])),col=tCols,lwd=3)

  arrows(0,0,MB_real[i,1],MB_real[i,2],col="black",lwd=4,length=0.15)
  arrows(0,0,MB_pred[i,1],MB_pred[i,2],col=rgb(0.5,0.5,0.5),lwd=4,length=0.15)

  text(-sign(mean(cumsum(Re(seM[,2]))))*0.5,seq(-0.5,0.67,length.out=ncol(ChIP)),colnames(ChIP)[fOrd],col=tCols,cex=1.3,font=2)
  text(-0.65,-0.68,i,cex=1.3)
}
plot(0,0,type="n",axes=FALSE,xlim=c(-1.5,1.5),ylim=c(-4,1))
arrows(1,0.5,-1,0.5,col="black",lwd=4,length=0.15)
arrows(1,-0.5,-1,-0.5,,col=rgb(0.5,0.5,0.5),lwd=4,length=0.15)
text(-0.6,0.8,"Measured",pos=4,cex=1.3,font=2)
text(-0.6,0.8-1,"Predicted",pos=4,cex=1.3,font=2)
#dev.off()

# select the remaining clusters to plot in the supplementary material
clusToPlot2 <- sort(setdiff(order(rowSums(abs(ChIP_cl_comp)),decreasing=TRUE),clusToPlot1))

#pdf("plots/predict_ATAC_by_currTFs_additive_TfClusters_Suppl.pdf",pointsize=18,width=20,height=24)
par(mfrow=c(10,8),mar=c(0.5,0.5,0.5,0.5))
for(i in clusToPlot2){
  plot(0,0,xlim=c(-0.6,0.6),ylim=c(-0.6,0.6),type="n",xlab="",ylab="",axes=FALSE)
  axis(1,labels = FALSE)
  axis(2,labels = FALSE)
  box()
  plotCircles(c(0.2,0.4,0.6),lty=2,col="gray")

  fOrd <- order(abs(ChIP_cl_comp[i,]))
  seM <- cbind(c(0,ChIP_cl_comp[i,fOrd][-length(fOrd)]),ChIP_cl_comp[i,fOrd])
  tCols <- tf_cols[fOrd]
  segments(cumsum(Re(seM[,1])),cumsum(Im(seM[,1])),cumsum(Re(seM[,2])),cumsum(Im(seM[,2])),col=tCols,lwd=3)

  arrows(0,0,MB_real[i,1],MB_real[i,2],col="black",lwd=4,length=0.15)
  arrows(0,0,MB_pred[i,1],MB_pred[i,2],col=rgb(0.5,0.5,0.5),lwd=4,length=0.15)

  text(-sign(mean(cumsum(Re(seM[,2]))))*0.4,seq(-0.4,0.57,length.out=ncol(ChIP)),colnames(ChIP)[fOrd],col=tCols,cex=1.3,font=2)
  text(-0.55,-0.58,i,cex=1.3)
}
plot(0,0,type="n",axes=FALSE,xlim=c(-1.5,1.5),ylim=c(-4,1))
arrows(1,0.5,-1,0.5,col="black",lwd=4,length=0.15)
arrows(1,-0.5,-1,-0.5,,col=rgb(0.5,0.5,0.5),lwd=4,length=0.15)
text(-0.6,0.8,"Measured",pos=4,cex=1.3,font=2)
text(-0.6,0.8-1,"Predicted",pos=4,cex=1.3,font=2)
#dev.off()



# convert TF activity vectors to lenght=1 and compensate the ChIP enrichment values (has no impact on model)
# permute the angles of the TFs, for each peak separately
fm_coeffs_cplx_U <- fm_coeffs_cplx/abs(fm_coeffs_cplx)
ChIP_U <- t(t(ChIP)*abs(fm_coeffs_cplx))

set.seed(13)
fmPredM_rnd <- t(sapply(1:nrow(ChIP),function(i){sample(fm_coeffs_cplx_U)})) * ChIP_U

PLD <- cbind(pathLength=rowSums(abs(fmPredM)),predLength=abs(rowSums(fmPredM)),predLengthRND=abs(rowSums(fmPredM_rnd)))

#pdf("plots/predict_ATAC_by_currTFs_TFClustes_pathLength_vs_sum.pdf",pointsize=18,width=7,height=7)
plot(PLD[,c(1,2)],pch='.',cex=3,ylim=c(0,1.3),xlab="Total path length",ylab="Length of prediction vector",col="darkgray")
lines(lowess(PLD[,1],PLD[,2],f=0.1), col="blue",lwd=2)
lines(lowess(PLD[,1],PLD[,3],f=0.1), col="red",lwd=2)
legend(x="topleft",lwd=2, legend=c("Actual data","Randomized control"),lty=1,col=c("blue","red"),bty="n",title="Loess fit")
#dev.off()

