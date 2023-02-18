rm(list = ls())

library(RColorBrewer)
library(wsyn)

if(! "rstudioapi" %in% installed.packages()[,"Package"]){
  install.packages("rstudioapi")
}
library("rstudioapi")
setwd(dirname(dirname(getSourceEditorContext()$path)))

##### Preliminary Data Processing/Management #####
ADChem_raw<-read.csv("./data/altm_data.csv", stringsAsFactors = F)


ADChem_Dat<-ADChem_raw
ADChem_Dat<-ADChem_Dat[order(ADChem_Dat$SiteID, ADChem_Dat$ConfDate),]
ADChem_Dat$month<-as.character(substr(ADChem_Dat$ConfDate,6,7))
ADChem_Dat$month_abs<-as.character(substr(ADChem_Dat$ConfDate,1,7))

## Do some easy subsetting first to reduce number of computations needed -------
ADChem_Dat<-ADChem_Dat[!is.na(ADChem_Dat$ConfDate),] #remove any rows without a date
ADChem_Dat<-ADChem_Dat[!is.na(ADChem_Dat$LabID),] #remove any rows without a Lab ID
ADChem_BLANK<-ADChem_Dat[ADChem_Dat$SiteID=="BLANK",] #save the blanks separately
ADChem_Dat<-ADChem_Dat[ADChem_Dat$SiteID!="BLANK",]
ADChem_Dat<-ADChem_Dat[ADChem_Dat$Station==1,] #Remove unwanted stations
ADChem_Dat<-ADChem_Dat[!is.na(ADChem_Dat$depth),]
ADChem_Dat<-ADChem_Dat[ADChem_Dat$depth<=0.5,] #Remove unwanted depths

any(duplicated(ADChem_BLANK$ConfDate))

ADSites<-sort(unique(ADChem_Dat$SiteID))
ADMonths<-sort(unique(ADChem_Dat$month))
ADMonths_abs<-sort(unique(ADChem_Dat$month_abs))
ADCats<-rownames(ADChem_Dat)

## Do blank correction --------------------------------------------------------
ADChem_Cor<-ADChem_Dat
ADChem_Cor$cor<-rep(0, nrow(ADChem_Cor))
cor.cols<-8:24 #save indices of columns to correct
for(ii in 1:nrow(ADChem_Cor)){
  same.date<-ADChem_BLANK$ConfDate==ADChem_Cor$ConfDate[ii]
  if(sum(same.date)>0){
    ADChem_Cor[ii, cor.cols]<- ADChem_Cor[ii, cor.cols] - colMeans(ADChem_BLANK[same.date, cor.cols], na.rm=T)
    ADChem_Cor$cor[ii]<-1
  }
}

## Adjust site IDs ------------------------------------------------------------
ADChem_Cor$SiteID[ADChem_Cor$SiteID=="040747A"]<-"40747"
ADChem_Cor$SiteID[ADChem_Cor$SiteID=="040750A"]<-"40750"

# unique(ADChem_Cor$SiteID)
# Reformatting for coordination's sake NEEDS TO BE DONE BEFORE TESTING CLUSTERING FUNCTION ####
ADSites_Corr<-sort(unique(as.character(ADChem_Cor$SiteID)))

locstrue<-read.csv("./data/lake_coords.csv",stringsAsFactors = F)
locstrue<-locstrue[order(locstrue$SiteID),]

limelocs<-c("30172", "40905", "60182", "40576")

coordsSet<-data.frame(NA, nrow = 53, ncol = 3)
colnames(coordsSet)<-c("Site","X","Y")

for(i in 1:53){
  for(j in 1:3){
    coordsSet[i,j]<-locstrue[i,j]
  }
}

coords<-coordsSet[,c(2,3)]
coords<-coords[-which(coordsSet[,1] %in% limelocs),]
#print(coords)

ADSites_Corr<-coordsSet[-which(coordsSet[,1] %in% limelocs),1]

## write function to create matrices and replace NAs with monthly medians
makeDataMatrix<-function(df, variable){
  
  subdat<-df[,c("ConfDate","SiteID","month","month_abs",variable)]
  
  outmat<-matrix(NA, nrow=length(unique(df$SiteID)), ncol=length(unique(df$ConfDate)))
  sites<-sort(unique(df$SiteID))
  times<-sort(unique(df$ConfDate))
  #organize into matrix
  for(site in 1:length(sites)){
    for(tt in 1:length(times)){
      outmat[site,tt]<-mean(subdat[,variable][subdat$SiteID==sites[site] & subdat$ConfDate==times[tt]], na.rm=T)
    }
  }
  rownames(outmat)<-sites
  colnames(outmat)<-times
  
  #replace NAs with monthly median
  months<-substr(colnames(outmat),6,7)
  for(site in 1:length(sites)){
    for(tt in 1:length(times)){
      if(is.na(outmat[site,tt])){
        outmat[site,tt]<-median(outmat[site, months==months[tt]], na.rm=T)
      }
    }
  }
  return(outmat)
}
# Editing ADChem_Cor to drop limed sites
ADChem_Cor<-ADChem_Cor[-which(ADChem_Cor$SiteID %in% limelocs),]
ADChem_Cor$DOC[ADChem_Cor$DOC < 0.5] <- NA
ADChem_Cor$DOC[ADChem_Cor$DOC > 25] <-NA
# Convention altered to include Mat so Fluorine could continue to use F
DOCMat<-makeDataMatrix(ADChem_Cor, variable="DOC")
SO4Mat<-makeDataMatrix(ADChem_Cor, variable="SO4")
NO3Mat<-makeDataMatrix(ADChem_Cor, variable="NO3")
ClMat<-makeDataMatrix(ADChem_Cor, variable="Cl")
FMat<-makeDataMatrix(ADChem_Cor, variable="F")
ANCMat<-makeDataMatrix(ADChem_Cor, variable="ANC")
DICMat<-makeDataMatrix(ADChem_Cor, variable="DIC")
SiO2Mat<-makeDataMatrix(ADChem_Cor, variable="SiO2")
CaMat<-makeDataMatrix(ADChem_Cor, variable="Ca")
MgMat<-makeDataMatrix(ADChem_Cor, variable="Mg")
NaMat<-makeDataMatrix(ADChem_Cor, variable="Na")
KMat<-makeDataMatrix(ADChem_Cor, variable="K")
NH4Mat<-makeDataMatrix(ADChem_Cor, variable="NH4")
AL_TDMat<-makeDataMatrix(ADChem_Cor, variable="AL_TD")
AL_TMMat<-makeDataMatrix(ADChem_Cor, variable="AL_TM")
AL_OMMat<-makeDataMatrix(ADChem_Cor, variable="AL_OM")
AL_IMMat<-makeDataMatrix(ADChem_Cor, variable="AL_IM")


## NADP #####
DepRaw<-read.csv("./data/acid_deposition.csv")

DepSO4Mat<-matrix(DepRaw$DepSO4, nrow=49, ncol=304, byrow=TRUE)
DepNMat<-matrix(DepRaw$DepN, nrow=49, ncol=304, byrow=TRUE)


## Solar radiation ##
RadRaw <- read.csv("./data/solar_radiation.csv")

PARMat <- matrix(RadRaw$par, nrow=49, ncol=237, byrow=TRUE)
UVMat <- matrix(RadRaw$uv, nrow=49, ncol=237,)


# NAO
#Reading in data
NAO<-read.csv("./data/nao.csv")
print(NAO[42*12+6,])
print(NAO[42*12+6+303,])
NAO<-NAO[(42*12+6):(42*12+6+303),]


NAOMat<-matrix(NAO$NAO, nrow=length(ADSites_Corr),ncol=nrow(NAO), byrow=TRUE)
rownames(NAOMat)<-ADSites_Corr
colnames(NAOMat)<-paste(NAO$yyyy,NAO$mm,sep = "x")
#NAOMatsub<-NAOMat[,(48*12+1):(48*12+237)]
#Subsetting matrices to PAR time range and duping ADKPARMat for consistent formatting



Ppt<-read.csv("./data/precip.csv")
PptMat<-as.matrix(Ppt[,-1])


# Cleaning
times<-sort(unique(ADChem_Cor$ConfDate))
cltimes<-seq(1:length(times))
clDOC<-cleandat(DOCMat,cltimes,clev=5)
clSO4<-cleandat(SO4Mat,cltimes,clev=5)
clNO3<-cleandat(NO3Mat,cltimes,clev=5)
clCl<-cleandat(ClMat,cltimes,clev=5)
clF<-cleandat(FMat,cltimes,clev=5)
clANC<-cleandat(ANCMat,cltimes,clev=5)
clDIC<-cleandat(DICMat,cltimes,clev=5)
clSiO2<-cleandat(SiO2Mat,cltimes,clev=5)
clCa<-cleandat(CaMat,cltimes,clev=5)
clMg<-cleandat(MgMat,cltimes,clev=5)
clNa<-cleandat(NaMat,cltimes,clev=5)
clK<-cleandat(KMat,cltimes,clev=5)
clNH4<-cleandat(NH4Mat,cltimes,clev=5)
clAL_TD<-cleandat(AL_TDMat,cltimes,clev=5)
clAL_TM<-cleandat(AL_TMMat,cltimes,clev=5)
clAL_IM<-cleandat(AL_IMMat,cltimes,clev=5)
clAL_OM<-cleandat(AL_OMMat,cltimes,clev=5)
clDeppH<-cleandat(DeppHMat,cltimes,clev=5)
clDepSO4<-cleandat(DepSO4Mat,cltimes,clev=5)
clDepN<-cleandat(DepNMat,cltimes,clev=5)
clPpt<-cleandat(PptMat,cltimes,clev = 5)
clNAO <- cleandat(NAOMat,cltimes, clev=5)

cltimes2 <- 1:ncol(PARMat)
DOCMatsub<-DOCMat[,-c(1:67)]
DepNMatsub<-DepNMat[,-c(1:67)]
DepSO4Matsub<-DepSO4Mat[,-c(1:67)]
ANCMatsub<-ANCMat[,-c(1:67)]
PptMatsub<-PptMat[,-c(1:67)]
NAOMatsub<-NAOMat[,-c(1:67)]

clDOCsub <- cleandat(DOCMatsub, cltimes2, clev=5)
clDepNsub <- cleandat(DepNMatsub, cltimes2, clev=5)
clDepSO4sub <- cleandat(DepSO4Matsub, cltimes2, clev=5)
clANCsub <- cleandat(ANCMatsub, cltimes2, clev=5)
clPptsub <- cleandat(PptMatsub, cltimes2, clev=5)
clNAOsub <- cleandat(NAOMatsub, cltimes2, clev=5)
clPAR <- cleandat(PARMat, cltimes2, clev=5)
clUV <- cleandat(UVMat, cltimes2, clev=5)




## Nice DOC WPMF figure function

plotwmf<-function(object,zlims=NULL,neat=TRUE,colorfill=NULL,colorbar=TRUE,title=NULL,filename=NA,xlocs=NULL,xlabs=NULL,...)
{
  wav<-Mod(get_values(object))
  times<-get_times(object)
  timescales<-get_timescales(object)
  
  if(is.null(zlims)){
    zlims<-range(wav,na.rm=T)
  }else
  {
    rg<-range(wav,na.rm=T)
    if (rg[1]<zlims[1] || rg[2]>zlims[2])
    {
      stop("Error in plotmag.tts: zlims must encompass the z axis range of what is being plotted")
    }
  }
  if(neat){
    inds<-which(!is.na(colMeans(wav,na.rm=T)))
    wav<-wav[,inds]
    timescales<-timescales[inds]
  }
  if(is.null(colorfill)){
    jetcolors <- c("#00007F", "blue", "#007FFF", "cyan", 
                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
    colorfill<-grDevices::colorRampPalette(jetcolors)
  }
  ylocs<-pretty(timescales,n=8)
  if(is.null(xlocs)){
    xlocs<-pretty(times,n=8)
  }
  
  if(is.null(xlabs)){
    xlabs<-xlocs
  }
  
  if (!is.na(filename))
  {
    grDevices::pdf(paste0(filename,".pdf"))
  }
  if (!colorbar)
  {
    graphics::image(x=times,y=log2(timescales),z=wav,xlab="Time",zlim=zlims,
                    ylab="Timescale (months)",axes=F,col=colorfill(100),main=title,...)
    graphics::axis(1,at = xlocs,labels=xlabs)
    #graphics::axis(2,at = log2(ylocs),labels = ylocs)
  }else
  {
    fields::image.plot(x=times,y=log2(timescales),z=wav,xlab="Time",zlim=zlims,
                       ylab="Timescale (months)",axes=F,col=colorfill(100),main=title,...)
    graphics::axis(1,at = xlocs,labels=xlabs)
    graphics::axis(2,at = log2(ylocs),labels = ylocs)
  }
  if (!is.na(filename))
  {
    grDevices::dev.off()
  }
}


DOCres <- wpmf(clDOC$cdat, cltimes, sigmethod="quick", nrand=10000)

xlocs <- seq(8, 304, by=24)
xlabs <- seq(1993, 2017, by=2)
ylocs <- c(0,4,12,24,48,96,144)

jetcolors <- c("#00007F", "blue", "#007FFF", "cyan", 
               "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
colorfill<-grDevices::colorRampPalette(jetcolors)


## Figure 2 -- WPMF


tiff("./Graphics/Fig2_timeseries_wpmf.tif", units="in", res=300, width=3.25, height=5)

layout(matrix(c(1:3,6:4), nrow=3, byrow=F), widths=c(0.9,0.1))
par(mar=c(3,3,1.1,0.5), mgp=c(1.8,0.5,0), tcl=-0.3)
plot(cltimes, rep(NA,length.out=length(times)), ylim=c(0,25), xaxt="n", 
     ylab=expression(DOC~(mg~L^-1)), xlab="Time")
axis(1, at=xlocs, labels=xlabs)
pal <- rainbow(nrow(DOCMat))
for(ii in 1:nrow(DOCMat)){
  lines(cltimes, DOCMat[ii,], col=pal[ii])
}
text(par("usr")[1]+0.05*(diff(par("usr")[1:2])), par("usr")[4]-0.075*(diff(par("usr")[3:4])), "a)")

plot(cltimes, colMeans(DOCMat, na.rm=T), type="l", xaxt="n",
     ylab=expression(DOC~(mg~L^-1)), xlab="Time")
axis(1, at=xlocs, labels=xlabs)
text(par("usr")[1]+0.05*(diff(par("usr")[1:2])), par("usr")[4]-0.075*(diff(par("usr")[3:4])), "b)")

plotwmf(DOCres, xlocs=xlocs, xlabs=xlabs, colorbar=FALSE)
axis(2, at=log2(ylocs), labels=ylocs)
q<-stats::quantile(DOCres$signif[[2]], 0.999)
graphics::contour(x=DOCres$times,y=log2(DOCres$timescales),z=Mod(DOCres$values),levels=q,drawlabels=F,lwd=1.5,
                  xaxs="i",xaxt="n",yaxt="n",xaxp=c(0,1,5),las = 1,frame=F, add=T)
text(par("usr")[1]+0.05*(diff(par("usr")[1:2])), par("usr")[4]-0.075*(diff(par("usr")[3:4])), "c)")

par(mar=c(3.6,1.1,1.6,0.5))
image(t(matrix(1:100)), col=colorfill(100), xaxt="n", yaxt="n")
axis(2, at=seq(0,1,length.out=5))

dev.off()



### Coherence testing #############################################################################

subann <- c(2,4)
annual <- c(9,16)
intann <- c(24,96)

getCohMag<-function(obj,band){
  x<-obj$coher[obj$timescales>=band[1] & obj$timescales<=band[2]]
  mnmag<-Mod(mean(x, na.rm=T))
  return(mnmag)
}

DOCxF <- coh(clDOC$cdat, clF$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxF <- bandtest(DOCxF, subann)
DOCxF <- bandtest(DOCxF, annual)
DOCxF <- bandtest(DOCxF, intann)
plotmag(DOCxF); mtext("DOCxF")
plotphase(DOCxF); mtext("DOCxF")
getCohMag(DOCxF, subann)
getCohMag(DOCxF, annual)

DOCxCl <- coh(clDOC$cdat, clCl$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxCl <- bandtest(DOCxCl, subann)
DOCxCl<- bandtest(DOCxCl, annual)
DOCxCl <- bandtest(DOCxCl, intann)
plotmag(DOCxCl); mtext("DOCxCl")
plotphase(DOCxCl); mtext("DOCxCl")
getCohMag(DOCxCl, annual)
getCohMag(DOCxCl, intann)

DOCxSO4 <- coh(clDOC$cdat, clSO4$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxSO4 <- bandtest(DOCxSO4, subann)
DOCxSO4<- bandtest(DOCxSO4, annual)
DOCxSO4 <- bandtest(DOCxSO4, intann)
plotmag(DOCxSO4); mtext("DOCxSO4")
plotphase(DOCxSO4); mtext("DOCxSO4")
getCohMag(DOCxSO4, subann)
getCohMag(DOCxSO4, annual)
getCohMag(DOCxSO4, intann)

DOCxNO3 <- coh(clDOC$cdat, clNO3$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxNO3 <- bandtest(DOCxNO3, subann)
DOCxNO3<- bandtest(DOCxNO3, annual)
DOCxNO3 <- bandtest(DOCxNO3, intann)
plotmag(DOCxNO3); mtext("DOCxNO3")
plotphase(DOCxNO3); mtext("DOCxNO3")
getCohMag(DOCxNO3, annual)
getCohMag(DOCxNO3, intann)

DOCxAL_TD <- coh(clDOC$cdat, clAL_TD$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxAL_TD <- bandtest(DOCxAL_TD, subann)
DOCxAL_TD <- bandtest(DOCxAL_TD, annual)
DOCxAL_TD <- bandtest(DOCxAL_TD, intann)
plotmag(DOCxAL_TD); mtext("DOCxAL_TD")
plotphase(DOCxAL_TD); mtext("DOCxAL_TD")
getCohMag(DOCxAL_TD, subann)
getCohMag(DOCxAL_TD, annual)
getCohMag(DOCxAL_TD, intann)

DOCxAL_IM <- coh(clDOC$cdat, clAL_IM$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxAL_IM <- bandtest(DOCxAL_IM, subann)
DOCxAL_IM <- bandtest(DOCxAL_IM, annual)
DOCxAL_IM <- bandtest(DOCxAL_IM, intann)
plotmag(DOCxAL_IM); mtext("DOCxAL_IM")
plotphase(DOCxAL_IM); mtext("DOCxAL_IM")
getCohMag(DOCxAL_IM, subann)
getCohMag(DOCxAL_IM, annual)

DOCxK <- coh(clDOC$cdat, clK$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxK <- bandtest(DOCxK, subann)
DOCxK <- bandtest(DOCxK, annual)
DOCxK <- bandtest(DOCxK, intann)
plotmag(DOCxK); mtext("DOCxK")
plotphase(DOCxK); mtext("DOCxK")
getCohMag(DOCxK, subann)

DOCxCa <- coh(clDOC$cdat, clCa$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxCa <- bandtest(DOCxCa, subann)
DOCxCa <- bandtest(DOCxCa, annual)
DOCxCa <- bandtest(DOCxCa, intann)
plotmag(DOCxCa); mtext("DOCxCa")
plotphase(DOCxCa); mtext("DOCxCa")
getCohMag(DOCxCa, subann)
getCohMag(DOCxCa, intann)

DOCxNH4 <- coh(clDOC$cdat, clNH4$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxNH4 <- bandtest(DOCxNH4, subann)
DOCxNH4 <- bandtest(DOCxNH4, annual)
DOCxNH4 <- bandtest(DOCxNH4, intann)
plotmag(DOCxNH4); mtext("DOCxNH4")
plotphase(DOCxNH4); mtext("DOCxNH4")
getCohMag(DOCxNH4, annual)

DOCxMg <- coh(clDOC$cdat, clMg$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxMg <- bandtest(DOCxMg, subann)
DOCxMg <- bandtest(DOCxMg, annual)
DOCxMg <- bandtest(DOCxMg, intann)
plotmag(DOCxMg); mtext("DOCxMg")
plotphase(DOCxMg); mtext("DOCxMg")
getCohMag(DOCxMg, subann)

DOCxNa <- coh(clDOC$cdat, clNa$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxNa <- bandtest(DOCxNa, subann)
DOCxNa <- bandtest(DOCxNa, annual)
DOCxNa <- bandtest(DOCxNa, intann)
plotmag(DOCxNa); mtext("DOCxNa")
plotphase(DOCxNa); mtext("DOCxNa")
getCohMag(DOCxNa, subann)

DOCxSiO2 <- coh(clDOC$cdat, clSiO2$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxSiO2 <- bandtest(DOCxSiO2, subann)
DOCxSiO2 <- bandtest(DOCxSiO2, annual)
DOCxSiO2 <- bandtest(DOCxSiO2, intann)
plotmag(DOCxSiO2); mtext("DOCxSiO2")
plotphase(DOCxSiO2); mtext("DOCxSiO2")
getCohMag(DOCxSiO2, annual)
getCohMag(DOCxSiO2, intann)

DOCxDIC <- coh(clDOC$cdat, clDIC$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxDIC <- bandtest(DOCxDIC, subann)
DOCxDIC <- bandtest(DOCxDIC, annual)
DOCxDIC <- bandtest(DOCxDIC, intann)
plotmag(DOCxDIC); mtext("DOCxDIC")
plotphase(DOCxDIC); mtext("DOCxDIC")
getCohMag(DOCxDIC, subann)
getCohMag(DOCxDIC, annual)

DOCxDepN <- coh(clDOC$cdat, clDepN$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxDepN <- bandtest(DOCxDepN, subann)
DOCxDepN <- bandtest(DOCxDepN, annual)
DOCxDepN <- bandtest(DOCxDepN, intann)
plotmag(DOCxDepN); mtext("DOCxDepN")
plotphase(DOCxDepN); mtext("DOCxDepN")
getCohMag(DOCxDepN, annual)

DOCxDepSO4 <- coh(clDOC$cdat, clDepSO4$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxDepSO4 <- bandtest(DOCxDepSO4, subann)
DOCxDepSO4 <- bandtest(DOCxDepSO4, annual)
DOCxDepSO4 <- bandtest(DOCxDepSO4, intann)
plotmag(DOCxDepSO4); mtext("DOCxDepSO4")
plotphase(DOCxDepSO4); mtext("DOCxDepSO4")
getCohMag(DOCxDepSO4, annual)

DOCxANC <- coh(clDOC$cdat, clANC$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxANC <- bandtest(DOCxANC, subann)
DOCxANC <- bandtest(DOCxANC, annual)
DOCxANC <- bandtest(DOCxANC, intann)
plotmag(DOCxANC); mtext("DOCxANC")
plotphase(DOCxANC); mtext("DOCxANC")
getCohMag(DOCxANC, subann)
getCohMag(DOCxANC, annual)

DOCxPpt <- coh(clDOC$cdat, clPpt$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxPpt <- bandtest(DOCxPpt, subann)
DOCxPpt <- bandtest(DOCxPpt, annual)
DOCxPpt <- bandtest(DOCxPpt, intann)
plotmag(DOCxPpt); mtext("DOCxPpt")
plotphase(DOCxPpt); mtext("DOCxPpt")
getCohMag(DOCxPpt, subann)
getCohMag(DOCxPpt, annual)

DOCxNAO <- coh(clDOC$cdat, clNAO$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DOCxNAO <- bandtest(DOCxNAO, subann)
DOCxNAO <- bandtest(DOCxNAO, annual)
DOCxNAO <- bandtest(DOCxNAO, intann)
plotmag(DOCxNAO); mtext("DOCxNAO")
plotphase(DOCxNAO); mtext("DOCxNAO")

DOCxPAR <- coh(clDOCsub$cdat, clPAR$cdat, times=cltimes2, norm="powall", sigmethod="fast", nrand=10000)
DOCxPAR <- bandtest(DOCxPAR, subann)
DOCxPAR <- bandtest(DOCxPAR, annual)
DOCxPAR <- bandtest(DOCxPAR, intann)
plotmag(DOCxPAR); mtext("DOCxPAR")
plotphase(DOCxPAR); mtext("DOCxPAR")
getCohMag(DOCxPAR, annual)

DOCxUV <- coh(clDOCsub$cdat, clUV$cdat, times=cltimes2, norm="powall", sigmethod="fast", nrand=10000)
DOCxUV <- bandtest(DOCxUV, subann)
DOCxUV <- bandtest(DOCxUV, annual)
DOCxUV <- bandtest(DOCxUV, intann)
plotmag(DOCxUV); mtext("DOCxUV")
plotphase(DOCxUV); mtext("DOCxUV")
getCohMag(DOCxUV, annual)


## coherences among external drivers
DepNxDepSO4 <- coh(clDepN$cdat, clDepSO4$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DepNxDepSO4 <- bandtest(DepNxDepSO4, subann)
DepNxDepSO4 <- bandtest(DepNxDepSO4, annual)
DepNxDepSO4 <- bandtest(DepNxDepSO4, intann)
plotmag(DepNxDepSO4); mtext("DepNxDepSO4")
plotphase(DepNxDepSO4); mtext("DepNxDepSO4")
getCohMag(DepNxDepSO4, subann)
getCohMag(DepNxDepSO4, annual)
getCohMag(DepNxDepSO4, intann)

DepNxANC <- coh(clDepN$cdat, clANC$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DepNxANC <- bandtest(DepNxANC, subann)
DepNxANC <- bandtest(DepNxANC, annual)
DepNxANC <- bandtest(DepNxANC, intann)
plotmag(DepNxANC); mtext("DepNxANC")
plotphase(DepNxANC); mtext("DepNxANC")
getCohMag(DepNxANC, subann)
getCohMag(DepNxANC, annual)
getCohMag(DepNxANC, intann)

DepNxPpt <- coh(clDepN$cdat, clPpt$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DepNxPpt <- bandtest(DepNxPpt, subann)
DepNxPpt <- bandtest(DepNxPpt, annual)
DepNxPpt <- bandtest(DepNxPpt, intann)
plotmag(DepNxPpt); mtext("DepNxPpt")
plotphase(DepNxPpt); mtext("DepNxPpt")
getCohMag(DepNxPpt, subann)
getCohMag(DepNxPpt, annual)
getCohMag(DepNxPpt, intann)

DepNxNAO <- coh(clDepN$cdat, clNAO$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DepNxNAO <- bandtest(DepNxNAO, subann)
DepNxNAO <- bandtest(DepNxNAO, annual)
DepNxNAO <- bandtest(DepNxNAO, intann)
plotmag(DepNxNAO); mtext("DepNxNAO")
plotphase(DepNxNAO); mtext("DepNxNAO")
getCohMag(DepNxNAO, subann)
getCohMag(DepNxNAO, annual)
getCohMag(DepNxNAO, intann)

DepNxUV <- coh(clDepNsub$cdat, clUV$cdat, times=cltimes2, norm="powall", sigmethod="fast", nrand=10000)
DepNxUV <- bandtest(DepNxUV, subann)
DepNxUV <- bandtest(DepNxUV, annual)
DepNxUV <- bandtest(DepNxUV, intann)
plotmag(DepNxUV); mtext("DepNxUV")
plotphase(DepNxUV); mtext("DepNxUV")
getCohMag(DepNxUV, subann)
getCohMag(DepNxUV, annual)
getCohMag(DepNxUV, intann)

DepNxPAR <- coh(clDepNsub$cdat, clPAR$cdat, times=cltimes2, norm="powall", sigmethod="fast", nrand=10000)
DepNxPAR <- bandtest(DepNxPAR, subann)
DepNxPAR <- bandtest(DepNxPAR, annual)
DepNxPAR <- bandtest(DepNxPAR, intann)
plotmag(DepNxPAR); mtext("DepNxPAR")
plotphase(DepNxPAR); mtext("DepNxPAR")
getCohMag(DepNxPAR, subann)
getCohMag(DepNxPAR, annual)
getCohMag(DepNxPAR, intann)

DepSO4xANC <- coh(clDepSO4$cdat, clANC$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DepSO4xANC <- bandtest(DepSO4xANC, subann)
DepSO4xANC <- bandtest(DepSO4xANC, annual)
DepSO4xANC <- bandtest(DepSO4xANC, intann)
plotmag(DepSO4xANC); mtext("DepSO4xANC")
plotphase(DepSO4xANC); mtext("DepSO4xANC")
getCohMag(DepSO4xANC, subann)
getCohMag(DepSO4xANC, annual)
getCohMag(DepSO4xANC, intann)

DepSO4xPpt <- coh(clDepSO4$cdat, clPpt$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DepSO4xPpt <- bandtest(DepSO4xPpt, subann)
DepSO4xPpt <- bandtest(DepSO4xPpt, annual)
DepSO4xPpt <- bandtest(DepSO4xPpt, intann)
plotmag(DepSO4xPpt); mtext("DepSO4xPpt")
plotphase(DepSO4xPpt); mtext("DepSO4xPpt")
getCohMag(DepSO4xPpt, subann)
getCohMag(DepSO4xPpt, annual)
getCohMag(DepSO4xPpt, intann)

DepSO4xNAO <- coh(clDepSO4$cdat, clNAO$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
DepSO4xNAO <- bandtest(DepSO4xNAO, subann)
DepSO4xNAO <- bandtest(DepSO4xNAO, annual)
DepSO4xNAO <- bandtest(DepSO4xNAO, intann)
plotmag(DepSO4xNAO); mtext("DepSO4xNAO")
plotphase(DepSO4xNAO); mtext("DepSO4xNAO")
getCohMag(DepSO4xNAO, subann)
getCohMag(DepSO4xNAO, annual)
getCohMag(DepSO4xNAO, intann)

DepSO4xUV <- coh(clDepSO4sub$cdat, clUV$cdat, times=cltimes2, norm="powall", sigmethod="fast", nrand=10000)
DepSO4xUV <- bandtest(DepSO4xUV, subann)
DepSO4xUV <- bandtest(DepSO4xUV, annual)
DepSO4xUV <- bandtest(DepSO4xUV, intann)
plotmag(DepSO4xUV); mtext("DepSO4xUV")
plotphase(DepSO4xUV); mtext("DepSO4xUV")
getCohMag(DepSO4xUV, subann)
getCohMag(DepSO4xUV, annual)
getCohMag(DepSO4xUV, intann)

DepSO4xPAR <- coh(clDepSO4sub$cdat, clPAR$cdat, times=cltimes2, norm="powall", sigmethod="fast", nrand=10000)
DepSO4xPAR <- bandtest(DepSO4xPAR, subann)
DepSO4xPAR <- bandtest(DepSO4xPAR, annual)
DepSO4xPAR <- bandtest(DepSO4xPAR, intann)
plotmag(DepSO4xPAR); mtext("DepSO4xPAR")
plotphase(DepSO4xPAR); mtext("DepSO4xPAR")
getCohMag(DepSO4xPAR, subann)
getCohMag(DepSO4xPAR, annual)
getCohMag(DepSO4xPAR, intann)

ANCxPpt <- coh(clANC$cdat, clPpt$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
ANCxPpt <- bandtest(ANCxPpt, subann)
ANCxPpt <- bandtest(ANCxPpt, annual)
ANCxPpt <- bandtest(ANCxPpt, intann)
plotmag(ANCxPpt); mtext("ANCxPpt")
plotphase(ANCxPpt); mtext("ANCxPpt")
getCohMag(ANCxPpt, subann)
getCohMag(ANCxPpt, annual)
getCohMag(ANCxPpt, intann)

ANCxNAO <- coh(clANC$cdat, clNAO$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
ANCxNAO <- bandtest(ANCxNAO, subann)
ANCxNAO <- bandtest(ANCxNAO, annual)
ANCxNAO <- bandtest(ANCxNAO, intann)
plotmag(ANCxNAO); mtext("ANCxNAO")
plotphase(ANCxNAO); mtext("ANCxNAO")
getCohMag(ANCxNAO, subann)
getCohMag(ANCxNAO, annual)
getCohMag(ANCxNAO, intann)

ANCxUV <- coh(clANCsub$cdat, clUV$cdat, times=cltimes2, norm="powall", sigmethod="fast", nrand=10000)
ANCxUV <- bandtest(ANCxUV, subann)
ANCxUV <- bandtest(ANCxUV, annual)
ANCxUV <- bandtest(ANCxUV, intann)
plotmag(ANCxUV); mtext("ANCxUV")
plotphase(ANCxUV); mtext("ANCxUV")
getCohMag(ANCxUV, subann)
getCohMag(ANCxUV, annual)
getCohMag(ANCxUV, intann)

ANCxPAR <- coh(clANCsub$cdat, clPAR$cdat, times=cltimes2, norm="powall", sigmethod="fast", nrand=10000)
ANCxPAR <- bandtest(ANCxPAR, subann)
ANCxPAR <- bandtest(ANCxPAR, annual)
ANCxPAR <- bandtest(ANCxPAR, intann)
plotmag(ANCxPAR); mtext("ANCxPAR")
plotphase(ANCxPAR); mtext("ANCxPAR")
getCohMag(ANCxPAR, subann)
getCohMag(ANCxPAR, annual)
getCohMag(ANCxPAR, intann)

PptxNAO <- coh(clPpt$cdat, clNAO$cdat, times=cltimes, norm="powall", sigmethod="fast", nrand=10000)
PptxNAO <- bandtest(PptxNAO, subann)
PptxNAO <- bandtest(PptxNAO, annual)
PptxNAO <- bandtest(PptxNAO, intann)
plotmag(PptxNAO); mtext("PptxNAO")
plotphase(PptxNAO); mtext("PptxNAO")
getCohMag(PptxNAO, subann)
getCohMag(PptxNAO, annual)
getCohMag(PptxNAO, intann)

PptxUV <- coh(clPptsub$cdat, clUV$cdat, times=cltimes2, norm="powall", sigmethod="fast", nrand=10000)
PptxUV <- bandtest(PptxUV, subann)
PptxUV <- bandtest(PptxUV, annual)
PptxUV <- bandtest(PptxUV, intann)
plotmag(PptxUV); mtext("PptxUV")
plotphase(PptxUV); mtext("PptxUV")
getCohMag(PptxUV, subann)
getCohMag(PptxUV, annual)
getCohMag(PptxUV, intann)

PptxPAR <- coh(clPptsub$cdat, clPAR$cdat, times=cltimes2, norm="powall", sigmethod="fast", nrand=10000)
PptxPAR <- bandtest(PptxPAR, subann)
PptxPAR <- bandtest(PptxPAR, annual)
PptxPAR <- bandtest(PptxPAR, intann)
plotmag(PptxPAR); mtext("PptxPAR")
plotphase(PptxPAR); mtext("PptxPAR")
getCohMag(PptxPAR, subann)
getCohMag(PptxPAR, annual)
getCohMag(PptxPAR, intann)

NAOxUV <- coh(clNAOsub$cdat, clUV$cdat, times=cltimes2, norm="powall", sigmethod="fast", nrand=10000)
NAOxUV <- bandtest(NAOxUV, subann)
NAOxUV <- bandtest(NAOxUV, annual)
NAOxUV <- bandtest(NAOxUV, intann)
plotmag(NAOxUV); mtext("NAOxUV")
plotphase(NAOxUV); mtext("NAOxUV")
getCohMag(NAOxUV, subann)
getCohMag(NAOxUV, annual)
getCohMag(NAOxUV, intann)

NAOxPAR <- coh(clNAOsub$cdat, clPAR$cdat, times=cltimes2, norm="powall", sigmethod="fast", nrand=10000)
NAOxPAR <- bandtest(NAOxPAR, subann)
NAOxPAR <- bandtest(NAOxPAR, annual)
NAOxPAR <- bandtest(NAOxPAR, intann)
plotmag(NAOxPAR); mtext("NAOxPAR")
plotphase(NAOxPAR); mtext("NAOxPAR")
getCohMag(NAOxPAR, subann)
getCohMag(NAOxPAR, annual)
getCohMag(NAOxPAR, intann)

UVxPAR <- coh(clUV$cdat, clPAR$cdat, times=cltimes2, norm="powall", sigmethod="fast", nrand=10000)
UVxPAR <- bandtest(UVxPAR, subann)
UVxPAR <- bandtest(UVxPAR, annual)
UVxPAR <- bandtest(UVxPAR, intann)
plotmag(UVxPAR); mtext("UVxPAR")
plotphase(UVxPAR); mtext("UVxPAR")
getCohMag(UVxPAR, subann)
getCohMag(UVxPAR, annual)
getCohMag(UVxPAR, intann)

### Wavelet linear models #########################################################################



# Short timescale

STSet<-list(DOC=clDOC$cdat, Precip=clPpt$cdat, ANC=clANC$cdat)

STwlm<-wlm(STSet,cltimes,resp=1,pred=2:3,norm="powall", scale.max.input = 96)

STse<-syncexpl(STwlm)
STseVals<-STse[STse$timescales>=2 & STse$timescales<=4,]
round(100*colMeans(STseVals[,c(3:9)])/mean(STseVals$sync),4)


# Medium timescale

LTSet<-list(DOC=clDOCsub$cdat, Precip=clPptsub$cdat, UV = clUV$cdat, PAR = clPAR$cdat, 
            NDep = clDepNsub$cdat, SO4Dep = clDepSO4sub$cdat)

LTwlm<-wlm(LTSet,cltimes2,resp=1,pred=c(2:3,5:6),norm="powall", scale.max.input = 96)

LTse<-syncexpl(LTwlm)
LTseVals<-LTse[LTse$timescales>=9 & LTse$timescales<=16,]
round(100*colMeans(LTseVals[,c(3:16)])/mean(LTseVals$sync),4)


LTwlm2<-wlm(LTSet,cltimes2,resp=1,pred=c(2,4:6),norm="powall", scale.max.input = 96)

LTse2<-syncexpl(LTwlm2)
LTseVals2<-LTse2[LTse2$timescales>=9 & LTse2$timescales<=16,]
round(100*colMeans(LTseVals2[,c(3:16)])/mean(LTseVals2$sync),4)


#stvals <- c(8.66, 7.22, 0.63, 0.81)
#stlabs <- c("Total","Ppt","ANC","Inter.")

ltvals <- c(79.46, 12.38, 50.49, 11.62, 19.35, -14.38)
ltlabs <- c("Total","Ppt","UV","Dep N","Dep SO4","Inter")


ltvals2 <- c(78.64, 12.54, 36.11, 11.42, 18.03, 0.53)
ltlabs2 <- c("Total","Ppt","PAR","Dep N","Dep SO4","Inter")

library(RColorBrewer)

pal <- brewer.pal(6,"Set2")

tiff("./Graphics/synchrExpl.tif", units="in", width=5.5, height=3.5, res=300)

par(mar=c(3,4,1.5,1.1))
bp<-barplot(ltvals, names.arg=NULL, ylab="Synchrony explained (%)", col=pal)
axis(1, at=bp[-5], labels=ltlabs[-5], col=NA)
axis(1, at=bp[5], labels=expression('Dep SO'[4]), col=NA, padj=0.2)

dev.off()


tiff("./Graphics/synchrExplAlt.tif", units="in", width=5.5, height=3.5, res=300)

par(mar=c(3,4,1.5,1.1))
bp<-barplot(ltvals2, names.arg=NULL, ylab="Synchrony explained (%)", col=pal)
axis(1, at=bp[-5], labels=ltlabs2[-5], col=NA)
axis(1, at=bp[5], labels=expression('Dep SO'[4]), col=NA, padj=0.2)

dev.off()





cor.test(ADChem_Cor$ANC, ADChem_Cor$Chlorophyll_a, use="pairwise.complete.obs")
