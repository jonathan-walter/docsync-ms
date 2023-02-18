### Analysis of distance decay/spatial structure in synchrony

rm(list = ls())

library(RColorBrewer)
library(ncf)
library(wsyn)
library(rgdal)

setwd(dirname(dirname(getSourceEditorContext()$path)))
##### Preliminary Data Processing/Management #####

ADChem_raw<-read.csv("./data/altm_data.csv",
                     stringsAsFactors = F)

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


years <- as.numeric(substr(ADMonths_abs,1,4))

cv <- function(x){
  return(sd(x, na.rm=T)/mean(x, na.rm=T))
}

DOCMat.sub <- DOCMat[, years %in% 1993:2016]

totalCV <- apply(DOCMat.sub, 1, cv)

years2 <- rep(1993:2016, each=12)
years3 <- 1993:2016
DOCMat.ann <- matrix(NA, nrow=nrow(DOCMat.sub), ncol=length(years3))
for(ii in 1:nrow(DOCMat.sub)){
  for(jj in 1:length(years3)){
    DOCMat.ann[ii,jj] <- median(DOCMat.sub[ii, years2==years3[jj]], na.rm=T)
  }
}


annCV <- apply(DOCMat.ann, 1, cv)

summary(annCV/totalCV)

tiff("./Graphics/CVs.tif", units="in", width=4, height=7, res=300)

par(mar=c(3,3,1.25,0.5), mgp=c(1.8,0.5,0), tcl=-0.3, mfrow=c(3,1))
hist(totalCV, main="", xlab="CV of all observations")
mtext("a)", at=par("usr")[1])
hist(annCV, main="", xlab="CV of annual medians")
mtext("b)", at=par("usr")[1])
hist(annCV/totalCV, main="", xlab="Annual CV:All CV")
mtext("c)", at=par("usr")[1])

dev.off()