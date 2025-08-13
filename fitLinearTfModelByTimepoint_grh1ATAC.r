# This R script contains the code to fit the time point specific linear transcription factor model to the GRH-1 acute depletion data.

library(splines)

# Smooth inferred TF activities using spline regression with half as many splines as time points
smoothCo <- function(Coe){
  nrSpl <- ceiling(ncol(Coe)/2) # use half as many splines than there are time points
  bsplineb <- bs(setNames(1:ncol(Coe),colnames(Coe)),nrSpl,intercept=TRUE)
  fmSpl <- lm(t(Coe) ~ 0 + bsplineb)
  t(bsplineb %*% fmSpl$coefficients)
}

# load master table with all the required information (GSE288867 on GEO:https://www.ncbi.nlm.nih.gov/geo/)
MT <- read.delim("data/GSE288867_ATAC_peaks_xy_linModelCont_allInfo.txt",check.names=FALSE,as.is=TRUE)
TF_enr <- as.matrix(MT[,grep("ChIP_",colnames(MT))])
colnames(TF_enr) <- gsub("ChIP_","",colnames(TF_enr))

# load GRH-1 acute depletion ATAC-seq timecourse data (GSE288865 on GEO:https://www.ncbi.nlm.nih.gov/geo/)
TL_grh1 <- as.matrix(read.delim("data/GSE288865_expr_peaks_norm.tab")[rownames(TF_enr),])
annot <- do.call(rbind,strsplit(colnames(TL_grh1),"_"))
rownames(annot) <- colnames(TL_grh1)
annot_times <- as.numeric(gsub("h","",annot[,3]))

# specify a string for each sample that denotes the condition (replace the string Pre in the annotation by EtOH)
condStr <- gsub("Pre","EtOH",paste(annot[,4],annot[,2],sep="_"))
condL <- split(1:nrow(annot),factor(condStr,levels=unique(condStr)))
# duplicate the time point zero from EtOH to also represent Aux replicate 1
condL[["1_Aux"]] <- c(1,condL[["1_Aux"]])

# calculate average the two replicate GRH-1 experiments
TLA_grh1 <- 0.5*cbind(TL_grh1[,condL[["2_EtOH"]]]+TL_grh1[,condL[["1_EtOH"]]],TL_grh1[,condL[["2_Aux"]]]+TL_grh1[,condL[["1_Aux"]]])
colnames(TLA_grh1) <- gsub("_2$","",colnames(TLA_grh1))

# mean normalze by one period in EtOH (this needs to be adjusted for every dataset)
TLAD_grh1 <- TLA_grh1 - rowMeans(TLA_grh1[,paste0("grh1_EtOH_",21:31,"h")])

# fit linear TF model to each sample individually
TF_coeffs <- lm(TLAD_grh1 ~ TF_enr + 0)$coefficients
rownames(TF_coeffs) <- colnames(TF_enr)

# specify the samples for all the conditions as well as time points
sa_etoh <- paste0("grh1_EtOH_",21:31,"h")
sa_aux <- gsub("EtOH","Aux",sa_etoh)
sa_time <- as.numeric(gsub("h","",do.call(rbind,strsplit(sa_etoh,"_"))[,3]))

# perform smoothing on the inferred transcription factor activities
TF_coeffs_etoh <- smoothCo(TF_coeffs[,sa_etoh])
TF_coeffs_aux <- smoothCo(TF_coeffs[,sa_aux])

#pdf("plots/grh1_linModelSpline_perTimePoint.pdf",width=16,height=7,pointsize=16)
par(mfrow=c(2,5),mar=c(4.5,4.5,1,1))
for(i in 1:ncol(TF_enr)){
  plot(sa_time,TF_coeffs_etoh[i,],type='l',lwd=2,col="blue",ylim=c(-0.25,0.25),xlab="Time (h)",ylab="TF activity (log2)")
  lines(sa_time,TF_coeffs_aux[i,],lwd=2,col="red")

  if(i == 5){legend(x="bottomleft",lwd=2, legend=c("grh-1 depletion EtOH","grh-1 depletion Auxin"),col=c("blue","red"),bty="n")}
  text(mean(sa_time),0.24,colnames(TF_enr)[i])
}
#dev.off()
